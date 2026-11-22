# ABOUTME: Tests of the well-posedness guard and admission: networks, faces, classes,
# ABOUTME: verdicts, messages, the oracles' agreement, and the macros, parser and copy rule.
using Random

@testset "the guard's faces follow the spec's list and the oracle's" begin
    for spec in MECHANISM_TEST_SPECS
        cm = spec.mechanism.catalytic
        faces = EnzymeRates.Derivation._faces(cm)
        @test faces[1] == 0
        @test [Set(EnzymeRates.Derivation._face_names(cm, f)) for f in faces] ==
              _testhelper_faces(cm)
        # No face zeroes a substrate with a product: the guard lists no mixed corner.
        @test !any(f -> (n = EnzymeRates.Derivation._face_names(cm, f);
                         any(in(cm.metabolites[cm.substrates]), n) &&
                         any(in(cm.metabolites[cm.products]), n)), faces)
    end
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
    @test length.(EnzymeRates.Derivation._faces.(
        (bibi.catalytic, inhibited.catalytic, terter.catalytic))) == (7, 14, 9)
    # P bound as a modifier is P taken up: it is zeroed with P, and P adds no modifier
    # face. Name `n` of the guard is metabolite `n` (S, then P).
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
    @test EnzymeRates.Derivation._faces(copy_P.catalytic) == UInt64[0, 2, 1]
    @test EnzymeRates.Derivation._face_names(copy_P.catalytic, UInt64(2)) == [:P]
end

@testset "a law the guard admits can be 0/0 at a mixed corner, scored as a 0 in a fit" begin
    # Census mechanism d7a09031d5dd2fc6 (ter-bi). At B = P = Q = 0 the rapid-equilibrium
    # segment {EBD, EBDP, EDP} has no form free of B and P, so the law is 0/0
    # there; the guard's faces never zero one substrate with every product (spec section
    # 8). At A = P = Q = 0, another such corner, the law is 0.
    am = @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q
        steps: begin
            (E + B ⇌ EB, EDP + B ⇌ EBDP, EP + B ⇌ EBP) :: OnlyA
            (E + P ⇌ EP, EB + P ⇌ EBP, EBD + P ⇌ EBDP,
             EQ + P ⇌ EPQ) :: EqualAI
            (E + Q <--> EQ, EP + Q <--> EPQ) :: EqualAI
            (EAB + D <--> EPQ, EB + D <--> EBD,
             EBP + D <--> EBDP, EP + D <--> EDP) :: OnlyA
            (EB + A <--> EAB, EBD + A <--> EPQ) :: OnlyA
            EQ + A ⇌ EAQ :: EqualAI
        end
    end
    @test EnzymeRates.Derivation._guard(am) ===
          EnzymeRates.Derivation._WELL_POSED
    params, Keq = _testhelper_random_independent_params(am; rng = Random.MersenneTwister(1))
    rate = compile_rate_equation(am)
    @test isnan(rate((A = 1.0, B = 0.0, D = 1.0, P = 0.0, Q = 0.0), params, Keq))
    @test rate((A = 0.0, B = 1.0, D = 1.0, P = 0.0, Q = 0.0), params, Keq) == 0.0
    # A fit scores the NaN as it scores the 0 a well-defined law predicts there: a zero rate
    # is refused, and a nonzero rate costs the same at the corner as at A = P = Q = 0, where
    # the law is 0.
    _testhelper_table(rate, A, B) =
        (group = [1, 1], Rate = [1.0, rate], A = [1.0, A], B = [1.0, B], D = [1.0, 1.0],
         P = [0.1, 0.0], Q = [0.1, 0.0])
    reference = _testhelper_row_reference(_testhelper_table(0.5, 1.0, 0.0))
    @test_throws(ErrorException("Zero rate at row 2: log(0) is undefined"),
                 FittingProblem(am, _testhelper_table(0.0, 1.0, 0.0), 1.0; reference))
    x = [log(params[p]) for p in parameters(am)]
    at_nan = FittingProblem(am, _testhelper_table(0.5, 1.0, 0.0), 1.0; reference)
    at_zero = FittingProblem(am, _testhelper_table(0.5, 0.0, 1.0), 1.0; reference)
    @test EnzymeRates.Fitting.loss!(x, at_nan) == EnzymeRates.Fitting.loss!(x, at_zero)
end

@testset "the law and ε networks: a node per form and copy, two transitions per edge" begin
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
    bind, iso = _testhelper_edge(r2, "E + S ⇌ ES"), _testhelper_edge(r2, "F <--> E")
    @test EnzymeRates.Derivation._excluded_steps(r2, 2)[bind] &&
          !EnzymeRates.Derivation._excluded_steps(r2, 1)[bind]
    @test EnzymeRates.Mechanisms.steps_in_active_conformation(r2) == trues(4)
    # The inactive law deletes both tagged steps, which strands every form but E.
    @test EnzymeRates.Mechanisms.steps_in_inactive_conformation(r2) == falses(4)
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
    # Names index the metabolites [S, P]: the binding takes up one S, its release gives
    # it off.
    active_uptake = _testhelper_transition(eps, bind, 1, true)
    active_release = _testhelper_transition(eps, bind, 1, false)
    @test eps.change[:, active_uptake] == [1, 0]
    @test (eps.takes[active_uptake], eps.takes[active_release]) == (UInt64(1), UInt64(0))
    mm = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> E + P
        end
    end
    law_mm, eps_mm = EnzymeRates.Derivation._law_network(mm),
        EnzymeRates.Derivation._eps_network(
            mm, EnzymeRates.Derivation._re_flags(mm.catalytic),
            EnzymeRates.Derivation._built_reading(mm))
    @test (law_mm.nforms, law_mm.ncopies, length(law_mm.from)) == (2, 1, 4)
    @test all(f -> getfield(law_mm, f) == getfield(eps_mm, f),
              fieldnames(EnzymeRates.Derivation._Network))
    s_bind = _testhelper_edge(mm, "E + S ⇌ ES")
    @test [law_mm.fast[t] for t in eachindex(law_mm.from) if law_mm.edge[t] == s_bind] ==
          [true, true]
    @test count(law_mm.fast) == 2
end

@testset "the guard's fast segments" begin
    # Each mechanism's law network, and its ε-network with every step rapid-equilibrium,
    # whose fast cycle turns the reaction over. The walk runs breadth-first along the fast
    # transitions from each present node no earlier walk reached; offsets count from a
    # segment's first node, and the cycle lists transitions.
    D = EnzymeRates.Derivation
    bibi = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EPQ <--> EQ + P
            EQ ⇌ E + Q
        end
    end
    # Forms E, EA, EAB, EPQ, EQ; offset rows A, B, P, Q.
    s = D._fast_segments(D._law_network(bibi))
    @test s.segment == [1, 1, 1, 2, 1]
    @test s.segments == [[1, 2, 5, 3], [4]]
    @test s.offset == [0 1 1 0 0; 0 0 1 0 0; 0 0 0 0 0; 0 0 0 0 1]
    @test s.cycle == Int[]
    s = D._fast_segments(D._eps_network(bibi, trues(5), D._built_reading(bibi)))
    @test s.segment == [1, 1, 1, 1, 1]
    @test s.segments == [[1, 2, 5, 3, 4]]
    @test s.offset == [0 1 1 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 1 1]
    @test s.cycle == [5, 1, 3, 9, 7]

    allosteric = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES      :: EqualAI
            ES <--> EP    :: EqualAI
            EP ⇌ E + P      :: EqualAI
            E + I ⇌ EI      :: OnlyA
        end
    end
    # Forms E, EI, EP, ES in each copy, the inactive one from node 5; offset rows S, P, I,
    # the metabolites in stored order. The law network holds no inactive EI (node 6),
    # which no walk reaches.
    s = D._fast_segments(D._law_network(allosteric))
    @test s.segment == [1, 1, 1, 1, 1, 0, 1, 1]
    @test s.segments == [[1, 2, 3, 4, 5, 7, 8]]
    @test s.offset == [0 0 0 1 0 0 0 1; 0 0 1 0 0 0 1 0; 0 1 0 0 0 0 0 0]
    @test s.cycle == Int[]
    s = D._fast_segments(
        D._eps_network(allosteric, trues(4), D._built_reading(allosteric)))
    @test s.segment == [1, 1, 1, 1, 1, 1, 1, 1]
    @test s.segments == [[1, 2, 3, 4, 5, 6, 7, 8]]
    @test s.offset == [0 0 0 1 0 0 0 1; 0 0 1 0 0 0 1 0; 0 1 0 0 0 1 0 0]
    @test s.cycle == [3, 5, 8]

    pingpong = @_testhelper_unguarded @enzyme_mechanism begin
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
    # Forms E, EA, EB_res_plusA_minusP, EP_res_plusA_minusP, EQ, E_res_plusA_minusP;
    # offset rows A, B, P, Q. The residual segment is rooted at its first form,
    # EB_res_plusA_minusP, so E_res_plusA_minusP's offset of B is -1. With every step fast
    # the walk reaches E_res_plusA_minusP last, by the release of P, so its offset of P is
    # -1.
    s = D._fast_segments(D._law_network(pingpong))
    @test s.segment == [1, 1, 2, 2, 1, 2]
    @test s.segments == [[1, 2, 5], [3, 6, 4]]
    @test s.offset == [0 1 0 0 0 0; 0 0 0 -1 0 -1; 0 0 0 1 0 0; 0 0 0 0 1 0]
    @test s.cycle == Int[]
    s = D._fast_segments(D._eps_network(pingpong, trues(6), D._built_reading(pingpong)))
    @test s.segment == [1, 1, 1, 1, 1, 1]
    @test s.segments == [[1, 2, 5, 4, 3, 6]]
    @test s.offset == [0 1 0 1 0 1; 0 0 0 0 0 0; 0 0 0 0 0 -1; 0 0 1 0 1 0]
    @test s.cycle == [8, 3, 1, 5, 12, 10]
end

@testset "class 1 rejects a mechanism whose cycles never turn the reaction over" begin
    no_turnover = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP <--> F + P
        end
    end
    @test !EnzymeRates.Mechanisms.has_turnover_cycle(no_turnover.catalytic, trues(3))
    v = EnzymeRates.Derivation._guard(no_turnover)
    @test v === EnzymeRates.Derivation._Verdict(
        EnzymeRates.Derivation._NO_CYCLE, EnzymeRates.Derivation._LIMIT_UNTESTED, 0, 0, 0)
    @test EnzymeRates.Derivation._guard_message(no_turnover) ==
          "Mechanism: no cycle of steps turns substrates into products, so the rate is " *
          "zero everywhere."
    # With P released to free E the guard admits the mechanism, which has no message.
    turnover = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP <--> E + P
        end
    end
    @test EnzymeRates.Derivation._guard(turnover) === EnzymeRates.Derivation._WELL_POSED
    @test_throws "holds no fault" EnzymeRates.Derivation._guard_message(turnover)
    # The oracle's steady state carries no flux.
    setprecision(BigFloat, 256) do
        rng = Random.MersenneTwister(1)
        point = _testhelper_consistent_point(no_turnover, rng)
        net = _testhelper_network(no_turnover, point, Dict(:S => big"0.7", :P => big"1.3");
                                  speeds = _testhelper_draw_speeds(no_turnover, rng))
        @test abs(_testhelper_steady_rate(net; Λ = big(10.0)^20).v) < big(10)^-60
    end
end

@testset "the guard rejects a cycle that does not balance a metabolite" begin
    D = EnzymeRates.Derivation
    unbalanced = D._Verdict(D._UNBALANCED_CYCLE, D._LIMIT_UNTESTED, 0, 0, 0)
    rule = " Round every cycle, the substrates go in and the products come out as often " *
           "as the reaction runs, and each modifier comes off as often as it binds. " *
           "Check whether one of these steps leaves out a metabolite, writes one on the " *
           "wrong side or names the wrong form."
    # Form names are opaque, so a step that leaves out a metabolite makes no step wrong,
    # only a cycle. ES <--> ESI forgets `+ I`: the cycle E → EI → ESI → ES → E takes up
    # I and gives it off nowhere.
    forgotten_i = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES
            ES <--> E + P
            E + I ⇌ EI
            EI + S <--> ESI
            ES <--> ESI
        end
    end
    msg = "Mechanism: the steps E + S ⇌ ES, E + I ⇌ EI, EI + S <--> ESI and " *
          "ESI <--> ES close a cycle that, net, takes up I and gives off nothing." * rule
    @test D._guard(forgotten_i) === unbalanced
    @test !D.check_valid_rate_equation(Bool, forgotten_i)
    # The net change over S, P and I, read in the direction that takes up the first
    # metabolite it changes.
    @test D._unbalanced_cycle(forgotten_i.catalytic).change == [0, 0, 1]
    @test D._guard_message(forgotten_i) == msg
    @test_throws ErrorException(msg) parameters(forgotten_i)
    @test_throws ErrorException(msg) @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES
            ES <--> E + P
            E + I ⇌ EI
            EI + S <--> ESI
            ES <--> ESI
        end
    end
    # The same steps in an allosteric mechanism: the cycle is the catalytic scheme's.
    forgotten_i_allosteric = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES          :: EqualAI
            ES <--> E + P       :: OnlyA
            E + I ⇌ EI          :: EqualAI
            EI + S <--> ESI     :: EqualAI
            ES <--> ESI         :: EqualAI
        end
    end
    @test D._guard(forgotten_i_allosteric) === unbalanced
    @test !D.check_valid_rate_equation(Bool, forgotten_i_allosteric)
    @test D._guard_message(forgotten_i_allosteric) == msg
    # `EP + P ⇌ E` writes P's release as its uptake: the cycle E → ES → EP → E takes up S
    # and P.
    forgotten_release = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP + P ⇌ E
        end
    end
    msg = "Mechanism: the steps E + S ⇌ ES, EP + P ⇌ E and ES <--> EP close a cycle " *
          "that, net, takes up S and P and gives off nothing." * rule
    @test D._guard(forgotten_release) === unbalanced
    @test !D.check_valid_rate_equation(Bool, forgotten_release)
    @test D._unbalanced_cycle(forgotten_release.catalytic).change == [1, 1]
    @test D._guard_message(forgotten_release) == msg
    @test_throws ErrorException(msg) parameters(forgotten_release)
    @test_throws ErrorException(msg) @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP + P ⇌ E
        end
    end
    # EB ⇌ EAB forgets `+ A`, the first substrate: the cycle E → EA → EAB → EB → E takes
    # up A, while B goes on and comes off again. Measured against the uptake of A, as the
    # turnover count is, B would seem the one out of balance; the net change names A.
    forgotten_first = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EPQ ⇌ EQ + P
            EQ ⇌ E + Q
            E + B ⇌ EB
            EB ⇌ EAB
        end
    end
    msg = "Mechanism: the steps E + B ⇌ EB, E + A ⇌ EA, EA + B ⇌ EAB and EAB ⇌ EB close " *
          "a cycle that, net, takes up A and gives off nothing." * rule
    @test D._guard(forgotten_first) === unbalanced
    @test D._unbalanced_cycle(forgotten_first.catalytic).change == [1, 0, 0, 0]
    @test D._guard_message(forgotten_first) == msg
    @test_throws ErrorException(msg) @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EPQ ⇌ EQ + P
            EQ ⇌ E + Q
            E + B ⇌ EB
            EB ⇌ EAB
        end
    end
    # E + S + S ⇌ ESS takes up two S for each P, where the reaction takes up one: the
    # message counts the copies.
    two_s = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S + S ⇌ ESS
            ESS <--> E + P
        end
    end
    msg = "Mechanism: the steps E + S + S ⇌ ESS and E + P <--> ESS close a cycle that, " *
          "net, takes up 2 S and gives off P." * rule
    @test D._guard(two_s) === unbalanced
    @test D._unbalanced_cycle(two_s.catalytic).change == [2, -1]
    @test D._guard_message(two_s) == msg
    # E + I ⇌ EI and EI <--> E: the cycle takes up I and returns to E without it.
    consumes_i = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES
            ES <--> E + P
            E + I ⇌ EI
            EI <--> E
        end
    end
    msg = "Mechanism: the steps E <--> EI and E + I ⇌ EI close a cycle that, net, takes " *
          "up I and gives off nothing." * rule
    @test D._guard(consumes_i) === unbalanced
    # Read the other way, the cycle would give off I: the change is read in the direction
    # that takes up the first metabolite it changes, here I.
    @test D._unbalanced_cycle(consumes_i.catalytic).change == [0, 0, 1]
    @test D._guard_message(consumes_i) == msg
    @test_throws ErrorException(msg) @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES
            ES <--> E + P
            E + I ⇌ EI
            EI <--> E
        end
    end
    # The forgotten `+ I` written all rapid equilibrium: the fast cycle E → EI → ESI → ES
    # → E takes up I net, which class 2 would report as a fast cycle that turns the
    # reaction over, advising to make one of its steps steady-state; that advice leads to
    # `forgotten_i`. The check runs first and names the missing metabolite instead.
    forgotten_i_re = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES
            ES <--> E + P
            E + I ⇌ EI
            EI + S ⇌ ESI
            ESI ⇌ ES
        end
    end
    msg = "Mechanism: the steps E + S ⇌ ES, E + I ⇌ EI, EI + S ⇌ ESI and ESI ⇌ ES " *
          "close a cycle that, net, takes up I and gives off nothing." * rule
    @test D._guard(forgotten_i_re) === unbalanced
    @test D._guard_message(forgotten_i_re) == msg
    @test_throws ErrorException(msg) @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES
            ES <--> E + P
            E + I ⇌ EI
            EI + S ⇌ ESI
            ESI ⇌ ES
        end
    end
    # An allosteric mechanism whose unbalanced cycle, E → EI → E, lies among steps the
    # active conformation excludes: the check reads every cycle of the catalytic scheme,
    # whatever the conformations run, so the guard rejects it.
    consumes_i_inactive = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES        :: EqualAI
            ES <--> E + P     :: EqualAI
            E + I ⇌ EI        :: OnlyI
            EI <--> E         :: OnlyI
        end
    end
    msg = "Mechanism: the steps E <--> EI and E + I ⇌ EI close a cycle that, net, takes " *
          "up I and gives off nothing." * rule
    @test D._guard(consumes_i_inactive) === unbalanced
    active = EnzymeRates.Mechanisms.steps_in_active_conformation(consumes_i_inactive)
    @test !any(active[D._unbalanced_cycle(consumes_i_inactive.catalytic).steps])
    @test !D.check_valid_rate_equation(Bool, consumes_i_inactive)
    @test D._guard_message(consumes_i_inactive) == msg
    @test_throws ErrorException(msg) @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES        :: EqualAI
            ES <--> E + P     :: EqualAI
            E + I ⇌ EI        :: OnlyI
            EI <--> E         :: OnlyI
        end
    end
    # A cycle that runs the reaction twice takes up two S and gives off two P: balanced.
    twice = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES + S ⇌ ESS
            ESS <--> E + P + P
        end
    end
    @test D._guard(twice) === D._WELL_POSED
    @test abs.(D._fundamental_cycles(twice.catalytic,
                                     D._flat_edges(twice.catalytic))[2]) == [2]
    # S bound twice, by two orders: the cycle E → ES → ESSmod → ESmod → E takes up S twice
    # and gives it off twice.
    bound_twice = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> E + P
            ES + S ⇌ ESSmod
            E + S ⇌ ESmod
            ESmod + S ⇌ ESSmod
        end
    end
    @test D._guard(bound_twice) === D._WELL_POSED
end

@testset "class 2 rejects a rapid-equilibrium cycle that turns the reaction over" begin
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
    v = EnzymeRates.Derivation._guard(beside)
    @test (v.fault, v.conformation) == (EnzymeRates.Derivation._RE_TURNOVER, 1)
    @test EnzymeRates.Derivation._guard_message(beside) ==
          "Mechanism: the rapid-equilibrium steps E + P ⇌ EP, E + S ⇌ ES and " *
          "ES ⇌ EP close a cycle that turns the reaction over, so the mechanism has " *
          "no finite rate. Make one step of the cycle steady-state."
    beside_allo = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES      :: EqualAI
            ES ⇌ EP       :: EqualAI
            E + P ⇌ EP      :: EqualAI
            ES <--> FP    :: EqualAI
            FP <--> E + P   :: EqualAI
        end
    end
    allo_message = EnzymeRates.Derivation._guard_message(beside_allo)
    @test startswith(allo_message,
                     "Mechanism: the rapid-equilibrium steps E + P ⇌ EP, ")
    flip2 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            E + P ⇌ EP
            E + Q ⇌ EQ
            EA + B ⇌ EPQ
            EP + Q ⇌ EPQ
            EQ + P <--> EPQ
        end
    end
    # test_rate_eq_derivation.jl "All-RE catalytic cycle raises": its two mechanisms.
    allre = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S1, S2
        products: P1, P2
        steps: begin
            E + S1 ⇌ ES1
            E + S2 ⇌ ES2
            ES1 + S2 <--> ES1S2
            ES2 + S1 ⇌ ES1S2
            ES1S2 ⇌ EP1P2
            EP1P2 ⇌ EP1 + P2
            EP1P2 <--> EP2 + P1
            EP1 ⇌ E + P1
            EP2 ⇌ E + P2
        end
    end
    sib = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: R
        steps: begin
            E + S ⇌ ES
            ES ⇌ EP
            E + P <--> EP
            ER + S ⇌ ERS
            ERS ⇌ EPR
            ER + P ⇌ EPR
            (E + R ⇌ ER, ES + R ⇌ ERS, EP + R ⇌ EPR)
        end
    end
    fix = " Make one step of the cycle steady-state."
    for (i, m) in enumerate((beside, flip2, allre, sib))
        net = EnzymeRates.Derivation._eps_network(
            m, EnzymeRates.Derivation._re_flags(m.catalytic),
            EnzymeRates.Derivation._built_reading(m))
        cycle = EnzymeRates.Derivation._re_turnover_transitions(net)
        @test !isempty(cycle) && all(t -> net.fast[t], cycle)
        v = EnzymeRates.Derivation._guard(m)
        @test v.fault === EnzymeRates.Derivation._RE_TURNOVER
        msg = EnzymeRates.Derivation._guard_message(m)
        @test all(
            t -> occursin(EnzymeRates.Mechanisms.step_text(m.catalytic, net.edge[t]), msg),
            cycle)
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
                E + S ⇌ ES
                ES <--> EP
                E + P ⇌ EP
                ES <--> FP
                FP <--> E + P
            end
        end),
        @_testhelper_unguarded(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ EA
                E + P ⇌ EP
                E + Q ⇌ EQ
                EA + B <--> EPQ
                EP + Q ⇌ EPQ
                EQ + P <--> EPQ
            end
        end),
        @_testhelper_unguarded(@enzyme_mechanism begin
            substrates: S1, S2
            products: P1, P2
            steps: begin
                E + S1 ⇌ ES1
                E + S2 ⇌ ES2
                ES1 + S2 <--> ES1S2
                ES2 + S1 ⇌ ES1S2
                ES1S2 <--> EP1P2
                EP1P2 ⇌ EP1 + P2
                EP1P2 <--> EP2 + P1
                EP1 ⇌ E + P1
                EP2 ⇌ E + P2
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
            E + S ⇌ ES
            ES <--> EP
            E + P <--> EP
            ER + S ⇌ ERS
            ERS <--> EPR
            ER + P ⇌ EPR
            (E + R ⇌ ER, ES + R ⇌ ERS, EP + R ⇌ EPR)
        end
    end
    @test isempty(EnzymeRates.Derivation._re_turnover_transitions(
        EnzymeRates.Derivation._eps_network(
            sib_fixed, EnzymeRates.Derivation._re_flags(sib_fixed.catalytic),
            EnzymeRates.Derivation._built_reading(sib_fixed))))
end

@testset "the fitted classes pass classes 1 and 2" begin
    michaelis_menten = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> E + P
        end
    end
    abortive_complex = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B <--> EAB
            EAB <--> EPQ
            EQ + P <--> EPQ
            E + Q <--> EQ
            EQ + A ⇌ EAQ
        end
    end
    product_activated = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + P ⇌ EP
            EP + S ⇌ EPS
            EPS <--> EPP
            EP + P ⇌ EPP
        end
    end
    ss_flux_free_dead_end = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES
            ES <--> E + P
            E + I <--> EI
            ES <--> FS
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
            E + A ⇌ EA
            EA <--> EP_res_plusA_minusP
            EP_res_plusA_minusP ⇌ E_res_plusA_minusP + P
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
            EB_res_plusA_minusP <--> EQ
            EQ ⇌ E + Q
        end
    end
    theorell_chance = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B <--> EQ + P
            EQ ⇌ E + Q
        end
    end
    merged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B <--> EPQ
            EPQ ⇌ EQ + P
            EQ ⇌ E + Q
        end
    end
    copies = @_testhelper_unguarded @enzyme_mechanism begin
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
    mechs = [[s.mechanism for s in MECHANISM_TEST_SPECS];
             [pingpong, theorell_chance, merged, copies]]
    for mech in mechs
        cm = mech.catalytic
        g = _testhelper_graph(mech)
        names = (substrates = g.substrates, products = g.products,
                 modifiers = g.modifiers,
                 regulators = unique!([l for site in mech.allosteric.regulatory_sites
                                       for l in site.ligands]))
        for (k, e) in enumerate(cm.steps)
            op, lhs, rhs = EnzymeRates.DSL._step_parts(
                Meta.parse(EnzymeRates.Mechanisms.step_text(cm, k)))
            forms = Symbol[]
            s = EnzymeRates.DSL._raw_step!(forms, op, lhs, rhs, e.kinetic_group, names)
            @test forms == cm.enzyme_forms[[e.from, e.to]] &&
                  s == EnzymeRates.Mechanisms.Step(1, 2, e.taken_up, e.given_off,
                                                   e.rapid_equilibrium, e.kinetic_group)
        end
    end
    line(m, text) = EnzymeRates.Mechanisms.step_text(m.catalytic, _testhelper_edge(m, text))
    @test line(pingpong, "E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP") ==
          "E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP"
    @test line(theorell_chance, "EA + B <--> EQ + P") == "EA + B <--> EQ + P"
    @test line(merged, "EA + B <--> EPQ") == "EA + B <--> EPQ"
    @test line(copies, "ES + P ⇌ EPmodS") == "ES + P ⇌ EPmodS"
    @test line(copies, "E + I ⇌ EI") == "E + I ⇌ EI"
end

@testset "class 3 accepts the fitted classes and lumpable splits on every face" begin
    michaelis_menten = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> E + P
        end
    end
    abortive_complex = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B <--> EAB
            EAB <--> EPQ
            EQ + P <--> EPQ
            E + Q <--> EQ
            EQ + A ⇌ EAQ
        end
    end
    product_activated = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + P ⇌ EP
            EP + S ⇌ EPS
            EPS <--> EPP
            EP + P ⇌ EPP
        end
    end
    ss_flux_free_dead_end = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES
            ES <--> E + P
            E + I <--> EI
            ES <--> FS
        end
    end
    # R2's active conformation: a second free form, and no trap (spec section 8).
    r2_active = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            F + P ⇌ EP
            F <--> E
        end
    end
    tied_isoforms = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            (E + S <--> ES, F + S <--> FS)
            (ES <--> E + P, FS <--> F + P)
            (E + P ⇌ GP, F + P ⇌ GP)
            E <--> F
        end
    end
    untied_isoforms = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> ES
            F + S <--> FS
            (ES <--> E + P, FS <--> F + P)
            (E + P ⇌ GP, F + P ⇌ GP)
            E <--> F
        end
    end
    copy_trap = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> ES
            ES <--> EP
            EP ⇌ E + P
            EP + P ⇌ EPPmod
            EPPmod <--> EPmod + P
            E + P <--> EPmod
        end
    end
    lumpable_172 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ)
            (EA + B ⇌ EPQ, EQ + B ⇌ EBQ)
            (E + Q ⇌ EQ, EB + Q ⇌ EBQ)
            (EA + P <--> EAP, EQ + P <--> EPQ)
        end
    end
    lumpable_532 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            (EB + A ⇌ EPQ, EQ + A ⇌ EAQ)
            (E + B <--> EB, EA + B <--> EPQ, EI + B <--> EBI)
            (E + I ⇌ EI, EB + I ⇌ EBI, EQ + I ⇌ EIQ)
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ, EI + Q ⇌ EIQ)
            EQ + P <--> EPQ
        end
    end
    deadend_split = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> E + P
            E + P ⇌ EP
            F + P ⇌ EP
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
            E + A ⇌ EA
            EA + B <--> EAB
            EAB <--> EPQ
            EP + Q ⇌ EPQ
            EQ + P ⇌ EPQ
            E + P <--> EP
            E + Q ⇌ EQ
        end
    end
    faults = _testhelper_face_faults(ldh)
    @test faults == [ok, trap, ok, ok, ok, ok, split]
    v = EnzymeRates.Derivation._guard(ldh)
    @test (v.fault, v.limit, v.face, v.conformation, v.edge) ==
          (trap, EnzymeRates.Derivation._LIMIT_SPEED, 2, 1,
           _testhelper_edge(ldh, "EP + Q ⇌ EPQ"))
    @test EnzymeRates.Derivation._guard_message(ldh) ==
          "Mechanism: at zero P and Q, the rapid-equilibrium step EP + Q ⇌ EPQ " *
          "fills the forms {EP}, which no rapid-equilibrium step running there " *
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
            E + A ⇌ EA
            EA + B <--> EAB
            EAB <--> EPQ
            EP + Q <--> EPQ
            EQ + P ⇌ EPQ
            E + P <--> EP
            E + Q ⇌ EQ
        end
    end
    @test EnzymeRates.Derivation._guard(ldh_fed) === EnzymeRates.Derivation._WELL_POSED
    @test all(==(:agree),
              _testhelper_respeed_verdicts(ldh_fed; rng = Random.MersenneTwister(2)))
    child9 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA
            E + P ⇌ EP
            E + Q <--> EQ
            EA + B ⇌ EAB
            EA + Q <--> EAQ
            EP + B ⇌ EBP
            EP + Q ⇌ EPQ
            EQ + A <--> EAQ
            EQ + P ⇌ EPQ
            EAB <--> EPQ
        end
    end
    faults9 = _testhelper_face_faults(child9)
    @test faults9 == [ok, trap, ok, ok, ok, split, ok]
    v9 = _testhelper_face_verdicts(child9)[2]
    @test v9.edge == _testhelper_edge(child9, "EQ + P ⇌ EPQ")
    @test occursin("the rapid-equilibrium step EQ + P ⇌ EPQ fills the forms",
                   EnzymeRates.Derivation._face_message(child9, v9))
    verdicts9 = _testhelper_respeed_verdicts(child9; rng = Random.MersenneTwister(3))
    @test verdicts9[[2, 6]] == [:pattern, :pattern]
    @test all(_testhelper_agrees.(faults9, verdicts9))
    child9_fed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA
            E + P ⇌ EP
            E + Q <--> EQ
            EA + B ⇌ EAB
            EA + Q <--> EAQ
            EP + B ⇌ EBP
            EP + Q ⇌ EPQ
            EQ + A <--> EAQ
            EQ + P <--> EPQ
            EAB <--> EPQ
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
            E + A <--> EA
            E + B <--> EB
            E + Q ⇌ EQ
            EA + B <--> EPQ
            EB + Q ⇌ EBQ
            EB + A <--> EPQ
            EQ + B ⇌ EBQ
            EQ + P <--> EPQ
        end
    end
    faults = _testhelper_face_faults(split_child1)
    @test faults == [ok, split, ok, ok, ok, ok, split]
    at_Q = _testhelper_face_verdicts(split_child1)[7]
    msg = EnzymeRates.Derivation._face_message(split_child1, at_Q)
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
            E + A ⇌ EA
            EA + B <--> EAB
            EAB <--> EPQ
            EP + Q ⇌ EPQ
            EQ + P ⇌ EPQ
            E + P <--> EP
            E + Q ⇌ EQ
        end
    end
    @test EnzymeRates.Derivation._face_message(ldh,
                                               _testhelper_face_verdicts(ldh)[7]) ==
          "Mechanism: at zero Q, the forms {E, EA} and {EP} of one rapid-equilibrium " *
          "segment are joined only by steps that take up a zeroed metabolite " *
          "(E + Q ⇌ EQ and EP + Q ⇌ EPQ), so the rate law cannot weight one " *
          "against the other."
end

@testset "a split across both conformations names each form's conformation" begin
    # Both halves of the ping-pong :NonequalAI run chemistry in both conformations on one
    # rapid-equilibrium segment, joined through the flip at free E. At zero ATP and F6P
    # the covalent forms of each conformation make a part of their own, joined to the
    # rest only by steps that take up a zeroed metabolite. The two parts hold the same
    # forms, one part in each conformation, so the message names each form's
    # conformation and none for the whole split.
    am = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: ATP, F6P
        products: ADP, F16BP
        oligomeric_state: 1
        steps: begin
            E + ATP ⇌ EATP :: EqualAI
            EATP <--> EF16BP_res_plusATP_minusF16BP :: NonequalAI
            E_res_plusATP_minusF16BP + F16BP ⇌ EF16BP_res_plusATP_minusF16BP:: EqualAI
            E_res_plusATP_minusF16BP + F6P ⇌ EF6P_res_plusATP_minusF16BP :: EqualAI
            EF6P_res_plusATP_minusF16BP ⇌ EADP :: NonequalAI
            E + ADP ⇌ EADP :: EqualAI
        end
    end
    v = EnzymeRates.Derivation._guard(am)
    @test (v.fault, v.face, v.conformation) == (EnzymeRates.Derivation._SPLIT, 3, 0)
    @test EnzymeRates.Derivation._guard_message(am) ==
          "Mechanism: at zero ATP and F6P, the forms " *
          "{EF16BP_res_plusATP_minusF16BP (active), " *
          "E_res_plusATP_minusF16BP (active)} and " *
          "{EF16BP_res_plusATP_minusF16BP (inactive), " *
          "E_res_plusATP_minusF16BP (inactive)} of one rapid-equilibrium segment " *
          "are joined only by steps that take up a zeroed metabolite " *
          "(E_res_plusATP_minusF16BP + F6P ⇌ EF6P_res_plusATP_minusF16BP), so " *
          "the rate law cannot weight one against the other."
    # The LDH-like split at zero Q, every group :EqualAI: one part holds free E and EA
    # of both conformations, joined through the flip, so the labels go on the forms.
    ldh = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        oligomeric_state: 1
        steps: begin
            E + A ⇌ EA :: EqualAI
            EA + B <--> EAB :: EqualAI
            EAB <--> EPQ :: EqualAI
            EP + Q ⇌ EPQ :: EqualAI
            EQ + P ⇌ EPQ :: EqualAI
            E + P <--> EP :: EqualAI
            E + Q ⇌ EQ :: EqualAI
        end
    end
    at_Q = _testhelper_face_verdicts(ldh)[7]
    @test (at_Q.fault, at_Q.conformation) == (EnzymeRates.Derivation._SPLIT, 0)
    @test EnzymeRates.Derivation._face_message(ldh, at_Q) ==
          "Mechanism: at zero Q, the forms {E (active), EA (active), " *
          "E (inactive), EA (inactive)} and {EP (active)} and {EP (inactive)} of " *
          "one rapid-equilibrium segment are joined only by steps that take up a zeroed " *
          "metabolite (E + Q ⇌ EQ and EP + Q ⇌ EPQ), so the rate law cannot " *
          "weight one against the other."
end

@testset "a split whose parts lie in the inactive conformation is placed there" begin
    # A ping-pong whose X- and Y-bound routes run only in the inactive conformation. Free
    # E's segment starts at active E, but at zero P and Q its lowest parts, the X- and
    # Y-bound covalent forms, lie in the inactive conformation alone.
    am = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: X, Y
        steps: begin
            E + A ⇌ EA :: EqualAI
            EA <--> EP_res_plusA_minusP :: EqualAI
            EP_res_plusA_minusP ⇌ E_res_plusA_minusP + P :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            EB_res_plusA_minusP <--> EQ :: EqualAI
            EQ ⇌ E + Q :: EqualAI
            EA + X ⇌ EAX :: OnlyI
            EAX ⇌ EPX_res_plusA_minusP :: EqualAI
            EPX_res_plusA_minusP ⇌ EX_res_plusA_minusP + P :: EqualAI
            EX_res_plusA_minusP + B ⇌ EBX_res_plusA_minusP :: EqualAI
            EBX_res_plusA_minusP <--> EQX :: EqualAI
            EQ + X ⇌ EQX :: OnlyI
            EA + Y ⇌ EAY :: OnlyI
            EAY ⇌ EPY_res_plusA_minusP :: EqualAI
            EPY_res_plusA_minusP ⇌ EY_res_plusA_minusP + P :: EqualAI
            EY_res_plusA_minusP + B ⇌ EBY_res_plusA_minusP :: EqualAI
            EBY_res_plusA_minusP <--> EQY :: EqualAI
            EQ + Y ⇌ EQY :: OnlyI
        end
    end
    v = EnzymeRates.Derivation._guard(am)
    @test (v.fault, v.face, v.conformation) == (EnzymeRates.Derivation._SPLIT, 2, 2)
    @test EnzymeRates.Derivation._guard_message(am) ==
          "Mechanism: at zero P and Q in the inactive conformation, the forms " *
          "{EBX_res_plusA_minusP, EX_res_plusA_minusP} and " *
          "{EBY_res_plusA_minusP, EY_res_plusA_minusP} of one rapid-equilibrium " *
          "segment are joined only by steps that take up a zeroed metabolite " *
          "(EX_res_plusA_minusP + P ⇌ EPX_res_plusA_minusP and " *
          "EY_res_plusA_minusP + P ⇌ EPY_res_plusA_minusP), so the rate law " *
          "cannot weight one against the other."
end

@testset "a bottomless segment across both conformations labels each form" begin
    # At zero P and X the rapid-equilibrium segment of free E holds active E and, through
    # the flip, inactive E and the X-bound forms only the inactive conformation reaches.
    tc_re = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: X
        steps: begin
            E + S <--> ES          :: EqualAI
            ES <--> EP           :: EqualAI
            EP <--> E + P          :: EqualAI
            E + X ⇌ EX             :: OnlyI
            EX + S ⇌ FX + P      :: EqualAI
            FX <--> EX           :: EqualAI
        end
    end
    v = EnzymeRates.Derivation._guard(tc_re)
    @test (v.fault, v.face, v.conformation) == (EnzymeRates.Derivation._BOTTOMLESS, 5, 0)
    @test EnzymeRates.Derivation._guard_message(tc_re) ==
          "Mechanism: at zero P and X, the rapid-equilibrium segment " *
          "{E (active), E (inactive), EX (inactive), FX (inactive)} has no form " *
          "whose weight keeps the segment's lowest power of each zeroed metabolite, so " *
          "the forms' shares depend on how those concentrations vanish and the rate is " *
          "0/0 there."
    # No fix is printed: the segment's one binding made steady-state leaves a law with no
    # value on that face.
    tc_ss = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: X
        steps: begin
            E + S <--> ES          :: EqualAI
            ES <--> EP           :: EqualAI
            EP <--> E + P          :: EqualAI
            E + X <--> EX          :: OnlyI
            EX + S ⇌ FX + P      :: EqualAI
            FX <--> EX           :: EqualAI
        end
    end
    w = EnzymeRates.Derivation._guard(tc_ss)
    @test (w.fault, w.limit, w.face) ==
          (EnzymeRates.Derivation._LAW_MULTI, EnzymeRates.Derivation._LIMIT_LAW_NAN, 5)
end

@testset "a trap across both conformations labels each form" begin
    # Every step but the binding of B at EP runs in the active conformation alone. At
    # zero P and Q the fast class of free E holds E of both conformations, joined by the
    # flip, and the step E + P ⇌ EP that joins it to the rest of its segment takes up P.
    am = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + B <--> EB                                         :: OnlyA
            (E + P ⇌ EP, EB + P ⇌ EBP, EQ + P ⇌ EPQ)  :: OnlyA
            (E + Q <--> EQ, EP + Q <--> EPQ)                :: OnlyA
            EB + A ⇌ EPQ                                      :: OnlyA
            EP + B ⇌ EBP                                      :: EqualAI
        end
    end
    # EB + A ⇌ EPQ takes up A into EPQ, which the inactive conformation cuts
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
    @test EnzymeRates.Derivation._guard_message(am) ==
          "Mechanism: at zero P and Q, the forms {E (active), E (inactive)} " *
          "hold enzyme, but the rapid-equilibrium step E + P ⇌ EP that joins them to " *
          "the rest of their segment takes up a zeroed metabolite, so the rate law " *
          "misweights them. Make that step's kinetic group steady-state."
    # The face condition alone refuses the interior for that bridge; its message names
    # both tagged steps the cheapest tree needs.
    D = EnzymeRates.Derivation
    eps = D._eps_network(am, D._re_flags(am.catalytic), D._built_reading(am))
    w = D._face_condition(eps, D._law_network(am).present, UInt64(0),
                          first(am.catalytic.substrates))
    @test (w.fault, w.conformation, w.edge) ==
          (D._BRIDGE, 2, _testhelper_edge(am, "EB + P ⇌ EBP"))
    @test D._admission_message(am, D._Verdict(w.fault, w.limit, 1, w.conformation, w.edge),
                               D._built_reading(am)) ==
          "Mechanism: the :OnlyA steps EB + P ⇌ EBP; EQ + P ⇌ EPQ " *
          "bridge forms of the inactive conformation that its other steps keep apart, so " *
          "admission cannot show that the limit does not depend on how fast those steps " *
          "vanish. Make the kinetic group of EB + P ⇌ EBP steady-state."
    faults = _testhelper_face_faults(am)
    verdicts = _testhelper_respeed_verdicts(am; rng = Random.MersenneTwister(33))
    @test verdicts[[2, 6]] == [:mismatch, :mismatch]
    @test all(_testhelper_agrees.(faults, verdicts))
    # The fix: the binding of P steady-state.
    fixed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + B <--> EB                                                 :: OnlyA
            (E + P <--> EP, EB + P <--> EBP, EQ + P <--> EPQ) :: OnlyA
            (E + Q <--> EQ, EP + Q <--> EPQ)                        :: OnlyA
            EB + A ⇌ EPQ                                              :: OnlyA
            EP + B ⇌ EBP                                              :: EqualAI
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
            E + A ⇌ EA
            EA <--> EPQR
            EQR + P ⇌ EPQR
            EPR + Q ⇌ EPQR
            EPQ + R ⇌ EPQR
            ER + Q ⇌ EQR
            EQ + R ⇌ EQR
            ER + P ⇌ EPR
            EP + R ⇌ EPR
            EQ + P ⇌ EPQ
            EP + Q ⇌ EPQ
            E + P <--> EP
            E + Q <--> EQ
            E + R <--> ER
        end
    end
    v = EnzymeRates.Derivation._guard(terbot)
    @test (v.fault, v.face) == (EnzymeRates.Derivation._BOTTOMLESS, 2)
    @test EnzymeRates.Derivation._guard_message(terbot) ==
          "Mechanism: at zero P, Q and R, the rapid-equilibrium segment " *
          "{EP, EPQ, EPR, EPQR, EQ, ER, EQR} has no form whose " *
          "weight keeps the segment's lowest power of each zeroed metabolite, so the " *
          "forms' shares depend on how those concentrations vanish and the rate is 0/0 " *
          "there."
    terbot_ss = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A
        products: P, Q, R
        steps: begin
            E + A ⇌ EA
            EA <--> EPQR
            EQR + P ⇌ EPQR
            EPR + Q ⇌ EPQR
            EPQ + R ⇌ EPQR
            ER + Q ⇌ EQR
            EQ + R ⇌ EQR
            ER + P ⇌ EPR
            EP + R ⇌ EPR
            EQ + P ⇌ EPQ
            EP + Q <--> EPQ
            E + P <--> EP
            E + Q <--> EQ
            E + R <--> ER
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
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EP + Q ⇌ EPQ
            EQ + P ⇌ EPQ
            E + P <--> EP
            E + Q <--> EQ
        end
    end
    faults = _testhelper_face_faults(entered)
    @test faults == [ok, bottomless, ok, ok, ok, ok, ok]
    v = EnzymeRates.Derivation._guard(entered)
    @test (v.fault, v.face) == (bottomless, 2)
    @test EnzymeRates.Derivation._guard_message(entered) ==
          "Mechanism: at zero P and Q, the rapid-equilibrium segment " *
          "{EP, EPQ, EQ} has no form whose weight keeps the segment's lowest " *
          "power of each zeroed metabolite, so the forms' shares depend on how those " *
          "concentrations vanish and the rate is 0/0 there. Make the kinetic group of " *
          "EP + Q ⇌ EPQ steady-state."
    verdicts = _testhelper_respeed_verdicts(entered; rng = Random.MersenneTwister(6),
                                            law = :limit)
    @test verdicts[2] === :pattern
    @test all(_testhelper_agrees.(faults, verdicts))
    # The printed fix passes the guard.
    entered_named = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EP + Q <--> EPQ
            EQ + P ⇌ EPQ
            E + P <--> EP
            E + Q <--> EQ
        end
    end
    @test EnzymeRates.Derivation._guard(entered_named) ===
          EnzymeRates.Derivation._WELL_POSED
    # So does the segment's other binding made steady-state.
    entered_fixed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EP + Q ⇌ EPQ
            EQ + P <--> EPQ
            E + P <--> EP
            E + Q <--> EQ
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
            E + A ⇌ EA
            EA + B <--> EAB
            EAB <--> E + P + Q
            E + P <--> EP
            EP + Q ⇌ EPQ
            EQ + P ⇌ EPQ
            E + Q <--> EQ
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
            E + I ⇌ FI
            (E + P <--> EP, EP + P <--> EPP)
            E + S <--> EP
            E + S ⇌ FS
            EPP + S <--> FPP + P
        end
    end
    # At P = 0 the dead end FPS is closed and unreachable from E.
    unweighted = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> E + P
            E + P <--> EP
            EP + P <--> EPP
            EPP + S <--> FPS + P
        end
    end
    weighted = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> E + P
            E + P <--> EP
            EP + S <--> FS + P
        end
    end
    detour = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: X
        steps: begin
            E + S ⇌ ES              :: EqualAI
            ES <--> E + P           :: EqualAI
            E + P ⇌ EP              :: EqualAI
            EP + P ⇌ EPP        :: EqualAI
            FP + P <--> EPP     :: EqualAI
            E + X ⇌ EX              :: OnlyA
            FP + X <--> FPX     :: OnlyA
            EX + S ⇌ ESX        :: EqualAI
            ESX <--> FPX      :: EqualAI
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
    # The step between the conformations E and F, which the forms' first letters name.
    label(f) = first(String(weighted.catalytic.enzyme_forms[f]))
    leaving = only(k for (k, e) in enumerate(weighted.catalytic.steps)
                   if label(e.from) != label(e.to))
    @test v.edge == leaving
    @test EnzymeRates.Derivation._face_message(weighted, v) ==
          "Mechanism: at zero P, the forms {FS} keep enzyme away from free E: every " *
          "step that leaves them takes up a zeroed metabolite " *
          "($(EnzymeRates.Mechanisms.step_text(weighted.catalytic, leaving))), so the " *
          "steady state is not unique."
    # An allosteric message names the conformation the closed set lies in.
    at_P, at_PX = _testhelper_face_verdicts(detour)[[2, 5]]
    @test EnzymeRates.Derivation._face_message(detour, at_P) ==
          "Mechanism: at zero P in the inactive conformation, the forms {FP} " *
          "keep enzyme away from free E: every step that leaves them takes up a zeroed " *
          "metabolite (FP + P <--> EPP), so the steady state is not unique."
    @test EnzymeRates.Derivation._face_message(detour, at_PX) ==
          "Mechanism: at zero P and X in the active conformation, the forms " *
          "{FP} keep enzyme away from free E: every step that leaves them takes up a " *
          "zeroed metabolite (FP + P <--> EPP and FP + X <--> FPX), so the " *
          "steady state is not unique."
end

@testset "the limit test decides what structure cannot" begin
    unlumped = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA
            EB + A <--> EPQ
            (EA + B ⇌ EPQ, EQ + B ⇌ EBQ)
            (E + Q ⇌ EQ, EB + Q ⇌ EBQ)
            (EA + P <--> EAP, EQ + P <--> EPQ)
        end
    end
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
    dead_end = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P, Q
        steps: begin
            E + S ⇌ ES
            ES <--> EPQ
            EQ + P ⇌ EPQ
            E + Q <--> EQ
            EP + Q <--> EPQ
        end
    end
    # At zero P the aggregated chain has two closed classes, and so does the law's segment
    # chain, so the face stops at `_LAW_MULTI`.
    law_multi = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            E + P ⇌ EP
            EP <--> F + P
            F + S ⇌ FS
            FS <--> FP
            F + P ⇌ FP
            FP <--> E + P
        end
    end
    two_transients = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B <--> EAB
            EAB <--> EPQ
            EPQ ⇌ GPQ
            EQ + P ⇌ GPQ
            EP + Q ⇌ GPQ
            E + P <--> EP
            E + Q ⇌ EQ
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
    message(m, f) = EnzymeRates.Derivation._face_message(m,
                                                         _testhelper_face_verdicts(m)[f])
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
    differ = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EP_res_plusA_minusP                          :: EqualAI
            E + Q <--> EQ                                         :: NonequalAI
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP   :: NonequalAI
            E_res_plusA_minusP + B ⇌ EQ                        :: OnlyA
            EP_res_plusA_minusP + B ⇌ EBP_res_plusA_minusP  :: NonequalAI
            EQ + P ⇌ EBP_res_plusA_minusP                    :: NonequalAI
        end
    end
    # The :OnlyA binding of B into EQ, which the inactive conformation's Q binding still
    # reaches, cannot be switched off at rapid equilibrium, so the guard rejects the
    # interior first; class 3 sees the law chain at zero P and Q.
    v = EnzymeRates.Derivation._guard(differ)
    k = _testhelper_edge(differ, "E_res_plusA_minusP + B ⇌ EQ")
    @test v == EnzymeRates.Derivation._Verdict(
        EnzymeRates.Derivation._FAST_STEP_PINNED, EnzymeRates.Derivation._LIMIT_UNTESTED,
        0, 2, k)
    @test EnzymeRates.Derivation._guard_message(differ) ==
          _testhelper_pinned_cause(differ, k, 2)
    law = _testhelper_face_verdicts(differ)[2]
    @test law == EnzymeRates.Derivation._Verdict(
        EnzymeRates.Derivation._LAW_MULTI, EnzymeRates.Derivation._LIMIT_DIFFER, 2, 0, 0)
    @test EnzymeRates.Derivation._face_message(differ, law) ==
          "Mechanism: at zero P and Q, the steady state's fast-equilibrium " *
          "limit differs from the rate law."
    @test _testhelper_respeed_verdicts(differ; rng = Random.MersenneTwister(1),
                                       law = :limit)[2] === :mismatch
    # With that binding steady-state it is stopped, the interior passes, and the guard
    # itself reaches the law chain at zero P and Q; both oracles agree on every face.
    differ_ss = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EP_res_plusA_minusP                          :: EqualAI
            E + Q <--> EQ                                         :: NonequalAI
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP   :: NonequalAI
            E_res_plusA_minusP + B <--> EQ                     :: OnlyA
            EP_res_plusA_minusP + B ⇌ EBP_res_plusA_minusP  :: NonequalAI
            EQ + P ⇌ EBP_res_plusA_minusP                    :: NonequalAI
        end
    end
    v = EnzymeRates.Derivation._guard(differ_ss)
    @test v == EnzymeRates.Derivation._Verdict(
        EnzymeRates.Derivation._LAW_MULTI, EnzymeRates.Derivation._LIMIT_DIFFER, 2, 0, 0)
    @test EnzymeRates.Derivation._guard_message(differ_ss) ==
          "Mechanism: at zero P and Q, the steady state's fast-equilibrium " *
          "limit differs from the rate law."
    @test all(_testhelper_admission_agreement(differ_ss; rng = Random.MersenneTwister(2)))
end

@testset "a bridge's message names the tagged steps it needs" begin
    # E + A ⇌ EP_res_plusA_minusP binds A into a form the inactive conformation cuts
    # off, so it loses its affinity and stays fast: the cheapest spanning in-tree of that
    # conformation's fed forms needs only the vanishing :OnlyA binding of A at EQ, and
    # the message names it alone (the trap test's mechanism names a bridge of two). At zero
    # P and Q the ε-path converges: the rejection is conservative, as the message allows.
    bridge2 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EP_res_plusA_minusP                          :: OnlyA
            E + Q <--> EQ                                         :: EqualAI
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP   :: OnlyA
            E_res_plusA_minusP + B <--> EQ                     :: NonequalAI
            EP_res_plusA_minusP + Q ⇌ EAQ                    :: EqualAI
            EQ + A ⇌ EAQ                                      :: OnlyA
        end
    end
    v = EnzymeRates.Derivation._guard(bridge2)
    @test (v.fault, v.face, v.conformation, v.edge) ==
          (EnzymeRates.Derivation._BRIDGE, 1, 2, _testhelper_edge(bridge2, "EQ + A ⇌ EAQ"))
    @test EnzymeRates.Derivation._guard_message(bridge2) ==
          "Mechanism: the :OnlyA step EQ + A ⇌ EAQ bridges forms of the " *
          "inactive conformation that its other steps keep apart, so admission cannot " *
          "show that the limit does not depend on how fast that step vanishes. Make the " *
          "kinetic group of EQ + A ⇌ EAQ steady-state."
    @test _testhelper_epath_verdict(bridge2; rng = Random.MersenneTwister(2),
                                    face = Set([:P, :Q])) === :converges
    # With the binding of A into EP_res_plusA_minusP steady-state, the guard's first
    # built reading still finds a bridge at zero P and Q; the search admits the mechanism
    # under another reading, and the oracles agree on every face.
    leading_ss = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EP_res_plusA_minusP                       :: OnlyA
            E + Q <--> EQ                                         :: EqualAI
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP   :: OnlyA
            E_res_plusA_minusP + B <--> EQ                     :: NonequalAI
            EP_res_plusA_minusP + Q ⇌ EAQ                    :: EqualAI
            EQ + A ⇌ EAQ                                      :: OnlyA
        end
    end
    w = EnzymeRates.Derivation._admission(
        leading_ss, EnzymeRates.Derivation._built_reading(leading_ss))
    @test (w.fault, w.face, w.conformation) == (EnzymeRates.Derivation._BRIDGE, 2, 2)
    @test EnzymeRates.Derivation._guard(leading_ss) == EnzymeRates.Derivation._WELL_POSED
    @test all(_testhelper_admission_agreement(leading_ss; rng = Random.MersenneTwister(3),
                                              corners = true))
    _testhelper_check_law_against_limit(leading_ss)
end

@testset "the modifier faces catch a steady-state inhibitor binding" begin
    # At I = 0 the rapid-equilibrium square through EI and EIS no longer joins E
    # and ES, whose ratio the steady-state binding of S then sets.
    mixed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S <--> ES
            ES <--> E + P
            (E + I ⇌ EI, ES + I ⇌ EIS)
            EI + S ⇌ EIS
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
    # At zero I and S the slow reverse catalysis fills ES, and the rapid-equilibrium
    # binding of I that joins it to free E's segment takes up I.
    v = _testhelper_face_verdicts(mixed)[6]
    @test (v.fault, v.limit, v.conformation, v.edge) ==
          (EnzymeRates.Derivation._TRAP, EnzymeRates.Derivation._LIMIT_DIFFER, 1,
           _testhelper_edge(mixed, "ES + I ⇌ EIS"))
    @test EnzymeRates.Derivation._face_message(mixed, v) ==
          "Mechanism: at zero I and S, the forms {ES} hold enzyme, but the " *
          "rapid-equilibrium step ES + I ⇌ EIS that joins them to the rest of " *
          "their segment takes up a zeroed metabolite, so the rate law misweights them. " *
          "Make that step's kinetic group steady-state."
    # The fix: the binding of I steady-state.
    mixed_fixed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S <--> ES
            ES <--> E + P
            (E + I <--> EI, ES + I <--> EIS)
            EI + S ⇌ EIS
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
            (E + A <--> EA, EB + A <--> EAB, EP + A <--> EAP)
            (EA + B ⇌ EAB, EQ + B ⇌ EBQ)
            (E + P <--> EP, EA + P <--> EAP, EQ + P <--> EPQ)
            (E + Q ⇌ EQ, EB + Q ⇌ EBQ, EP + Q ⇌ EPQ)
            EAB <--> EPQ
        end
    end
    # The same mechanism with its groups and steps in reverse order, each binding
    # written as its release.
    rewritten = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            EAB <--> EPQ
            (EPQ ⇌ EP + Q, EBQ ⇌ EB + Q, EQ ⇌ E + Q)
            (EPQ <--> EQ + P, EAP <--> EA + P, EP <--> E + P)
            (EBQ ⇌ EQ + B, EAB ⇌ EA + B)
            (EAP <--> EP + A, EAB <--> EB + A, EA <--> E + A)
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
            E + A ⇌ EA
            EA + B <--> EAB
            EAB <--> EPQ
            EP + Q ⇌ EPQ
            EQ + P ⇌ EPQ
            E + P <--> EP
            E + Q ⇌ EQ
        end
    end
    candidate(re) = EnzymeRates.Enumeration._with_flags(_testhelper_candidate(ldh),
                                                        re).mechanism
    flags = _testhelper_flags(ldh)
    @test EnzymeRates.Derivation._guard(candidate(flags)).fault ===
          EnzymeRates.Derivation._TRAP
    # The feeding step steady-state through the flags: the graph of the trap test's fix.
    fed = copy(flags)
    fed[_testhelper_edge(ldh, "EP + Q ⇌ EPQ")] = false
    @test _testhelper_is_valid(candidate(fed))
    # Every step at rapid equilibrium: a fast cycle turns the reaction over.
    @test EnzymeRates.Derivation._guard(
              candidate(trues(length(ldh.catalytic.steps)))).fault ===
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
covalent form `E_res_plusA_minusP`, which gives the graph a free form `F` would;
AD2a_iso_chemNonequal is R1, AD2e is R2 and AD2c is R3 written again, and AD6 holds one
reaction in two groups, which the constructor rejects.
"""
const _testhelper_tag_sets = Dict{String, EnzymeRates.Mechanism}()

_testhelper_tag_sets["R1"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ ES        :: EqualAI
        ES <--> EP      :: NonequalAI
        F + P ⇌ EP        :: EqualAI
        F <--> E            :: NonequalAI
    end
end
_testhelper_tag_sets["A5"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ ES        :: EqualAI
        ES <--> EP      :: NonequalAI
        E + P ⇌ EP        :: EqualAI
        F ⇌ E               :: NonequalAI
    end
end
_testhelper_tag_sets["AD2a"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ ES        :: EqualAI
        ES <--> EP      :: OnlyA
        F + P ⇌ EP        :: EqualAI
        F <--> E            :: EqualAI
    end
end
_testhelper_tag_sets["AD2d"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ ES        :: EqualAI
        ES <--> EP      :: EqualAI
        F + P ⇌ EP        :: EqualAI
        F <--> E            :: OnlyA
    end
end
_testhelper_tag_sets["AD3a"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S ⇌ ES        :: EqualAI
        ES <--> EP      :: NonequalAI
        E + P ⇌ EP        :: EqualAI
        E + I ⇌ EI        :: OnlyA
    end
end
_testhelper_tag_sets["AD3b"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ ES                                       :: EqualAI
        ES <--> EP                                     :: EqualAI
        E + P ⇌ EP                                       :: EqualAI
        (E + S ⇌ ESmod, ES + S ⇌ ESSmod) :: OnlyA
        ESmod + S ⇌ ESSmod                       :: EqualAI
    end
end
_testhelper_tag_sets["AD3d"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ ES                    :: EqualAI
        ES <--> EP                  :: EqualAI
        E + P ⇌ EP                    :: EqualAI
        E + S <--> ESmod       :: OnlyA
        ES + S <--> ESSmod :: OnlyA
        ESmod + S ⇌ ESSmod    :: EqualAI
    end
end
_testhelper_tag_sets["AD4"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P, Q
    steps: begin
        E + S ⇌ ES            :: EqualAI
        ES <--> EPQ       :: EqualAI
        EPQ <--> E + P + Q  :: OnlyA
    end
end
_testhelper_tag_sets["AD4b"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P, Q
    steps: begin
        E + S ⇌ ES            :: EqualAI
        ES <--> EPQ       :: OnlyA
        EPQ ⇌ E + P + Q     :: EqualAI
    end
end
_testhelper_tag_sets["AD5"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    steps: begin
        E <--> F                        :: EqualAI
        (E + S ⇌ ES, F + S ⇌ FS)    :: OnlyA
        ES <--> FS                  :: EqualAI
        FS <--> FP                  :: OnlyA
        F + P ⇌ FP                    :: EqualAI
    end
end
_testhelper_tag_sets["AD7a"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: A
    products: P
    modifiers: I
    steps: begin
        E + A <--> EA                 :: EqualAI
        (E + I ⇌ EI, F + I ⇌ EI)    :: OnlyA
        E + P <--> EP                 :: EqualAI
        E ⇌ F                           :: EqualAI
        EA <--> EP                  :: EqualAI
    end
end
_testhelper_tag_sets["AD7b"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: A
    products: P
    modifiers: I
    steps: begin
        E + A <--> EA                 :: EqualAI
        (E + I ⇌ EI, F + I ⇌ EI)    :: EqualAI
        E + P <--> EP                 :: EqualAI
        E ⇌ F                           :: OnlyA
        EA <--> EP                  :: EqualAI
    end
end
_testhelper_tag_sets["AD7c"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: A
    products: P
    modifiers: I
    steps: begin
        E + A <--> EA                 :: EqualAI
        (E + I ⇌ EI, F + I ⇌ EI)    :: EqualAI
        E + P <--> EP                 :: EqualAI
        E ⇌ F                           :: EqualAI
        EA <--> EP                  :: OnlyA
    end
end
_testhelper_tag_sets["AD8"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A ⇌ EA            :: EqualAI
        EA + B <--> EQ + P  :: OnlyA
        E + Q ⇌ EQ            :: EqualAI
        EQ + B ⇌ EBQ      :: EqualAI
    end
end
_testhelper_tag_sets["AD9"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A ⇌ EA            :: EqualAI
        EA + B ⇌ EAB      :: EqualAI
        EAB ⇌ EPQ       :: OnlyA
        EQ + P <--> EPQ   :: EqualAI
        E + Q <--> EQ         :: EqualAI
        EA + Q <--> EAQ   :: OnlyA
    end
end
_testhelper_tag_sets["AD10"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A <--> EA         :: OnlyA
        E + B ⇌ EB            :: EqualAI
        EA + B ⇌ EAB      :: EqualAI
        EB + A <--> EAB   :: OnlyA
        EAB <--> EPQ    :: OnlyA
        EQ + P ⇌ EPQ      :: EqualAI
        E + Q ⇌ EQ            :: EqualAI
    end
end
# AD13: a second free form joined to E at rapid equilibrium, each label sorting before or
# after E, with the chemistry :NonequalAI or :OnlyA.
for label in (:D, :Ca, :Ea, :F, :Z), chem in (:NonequalAI, :OnlyA)
    block = quote
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES        :: EqualAI
            ES <--> EP      :: $chem
            E + P ⇌ EP        :: EqualAI
            $label ⇌ E          :: NonequalAI
        end
    end
    _testhelper_tag_sets["AD13_$(chem)_$label"] =
        EnzymeRates.DSL._construct(EnzymeRates.DSL.raw_mechanism(block))
end
_testhelper_tag_sets["AD14"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S ⇌ ES                            :: EqualAI
        E + P ⇌ EP                            :: EqualAI
        (ES <--> EP, EIS <--> EIP)  :: OnlyA
        EI + S ⇌ EIS                      :: OnlyA
        EI + P ⇌ EIP                      :: EqualAI
        E + I ⇌ EI                            :: EqualAI
    end
end
_testhelper_tag_sets["AD14b"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S ⇌ ES                            :: EqualAI
        E + P ⇌ EP                            :: EqualAI
        (ES <--> EP, EIS <--> EIP)  :: OnlyA
        EI + S <--> EIS                   :: OnlyA
        EI + P ⇌ EIP                      :: EqualAI
        E + I ⇌ EI                            :: EqualAI
    end
end
_testhelper_tag_sets["AD14c"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S ⇌ ES                            :: EqualAI
        E + P ⇌ EP                            :: EqualAI
        (ES <--> EP, EIS <--> EIP)  :: OnlyA
        EI + S ⇌ EIS                      :: EqualAI
        EI + P ⇌ EIP                      :: EqualAI
        E + I ⇌ EI                            :: EqualAI
    end
end
_testhelper_tag_sets["AD15"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        (E + A ⇌ EA, EB + A ⇌ EAB)  :: OnlyA
        EA + B ⇌ EAB                  :: OnlyA
        E + B ⇌ EB                        :: EqualAI
        EAB <--> EPQ                :: OnlyA
        EQ + P ⇌ EPQ                  :: EqualAI
        E + Q ⇌ EQ                        :: EqualAI
    end
end
_testhelper_tag_sets["R2P_sealed"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ ES        :: EqualAI
        ES <--> E + P     :: EqualAI
        E + P ⇌ EP        :: OnlyA
        F + P ⇌ EP        :: EqualAI
        F <--> E            :: OnlyA
    end
end
_testhelper_tag_sets["B_uphill"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    modifiers: X
    steps: begin
        E + S ⇌ ES        :: EqualAI
        ES <--> EP      :: EqualAI
        E + P ⇌ EP        :: EqualAI
        E + X ⇌ EX        :: OnlyA
        F + X <--> EX     :: EqualAI
        F <--> E            :: OnlyA
    end
end
_testhelper_tag_sets["TC1"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A ⇌ EA            :: EqualAI
        EA + B ⇌ EQ + P     :: OnlyA
        EQ <--> E + Q         :: EqualAI
        EA + B <--> EAB   :: EqualAI
        EAB <--> EPQ    :: EqualAI
        EPQ ⇌ EQ + P      :: EqualAI
    end
end
_testhelper_tag_sets["pp_V_half1"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A <--> EA                                       :: EqualAI
        EA <--> E_res_plusA_minusP + P                   :: OnlyA
        E_res_plusA_minusP + B <--> EB_res_plusA_minusP :: EqualAI
        EB_res_plusA_minusP <--> E + Q                     :: EqualAI
    end
end
_testhelper_tag_sets["pp_V_half2"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A <--> EA                                       :: EqualAI
        EA <--> E_res_plusA_minusP + P                   :: EqualAI
        E_res_plusA_minusP + B <--> EB_res_plusA_minusP :: EqualAI
        EB_res_plusA_minusP <--> E + Q                     :: OnlyA
    end
end
_testhelper_tag_sets["pp_K_Abind_half1"] =
    @_testhelper_unguarded @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A <--> EA                                       :: OnlyA
        EA <--> E_res_plusA_minusP + P                   :: OnlyA
        E_res_plusA_minusP + B <--> EB_res_plusA_minusP :: EqualAI
        EB_res_plusA_minusP <--> E + Q                     :: EqualAI
    end
end
_testhelper_tag_sets["PP1"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A ⇌ EA                                          :: EqualAI
        EA <--> E_res_plusA_minusP + P                   :: OnlyA
        E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP    :: EqualAI
        EB_res_plusA_minusP <--> E + Q                     :: EqualAI
    end
end
_testhelper_tag_sets["PP2"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    modifiers: I
    steps: begin
        E + A <--> EA                                       :: EqualAI
        EA <--> E_res_plusA_minusP + P                   :: EqualAI
        E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP    :: EqualAI
        EB_res_plusA_minusP <--> E + Q                     :: OnlyA
        E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP    :: OnlyA
    end
end
# Census ter-bi d75d2415361c4dae: its faces with P zero pass by T0+ alone.
_testhelper_tag_sets["terbi_d75d2415361c4dae"] =
    @_testhelper_unguarded @enzyme_mechanism begin
    substrates: A, B, D
    products: P, Q
    steps: begin
        (E + A ⇌ EA, EBP + A ⇌ EABP, EP + A ⇌ EAP)       :: OnlyA
        (E + P ⇌ EP, EA + P ⇌ EAP, EAB + P ⇌ EABP,
         EQ + P ⇌ EPQ)                                               :: EqualAI
        (E + Q <--> EQ, EP + Q <--> EPQ)                           :: OnlyA
        (EA + B <--> EAB, EAP + B <--> EABP,
         EP + B <--> EBP)                                            :: EqualAI
        (EAB + D <--> EPQ, EQ + D <--> EDQ)                  :: OnlyA
    end
end
# An :OnlyA group mixing B's fused binding into EPQ with its plain binding at EQ,
# which no untagged step joins to EBQ: its limit exists (spec section 7, "Mixed groups
# stay allowed").
_testhelper_tag_sets["mixed_group_limit"] =
    @_testhelper_unguarded @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A ⇌ EA                                    :: EqualAI
        (EA + B <--> EPQ, EQ + B <--> EBQ)  :: OnlyA
        EQ + P ⇌ EPQ                              :: EqualAI
        E + Q ⇌ EQ                                    :: EqualAI
    end
end
_testhelper_tag_sets["R2"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ ES        :: OnlyA
        ES <--> EP      :: EqualAI
        F + P ⇌ EP        :: EqualAI
        F <--> E            :: OnlyA
    end
end
_testhelper_tag_sets["BS1"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S ⇌ ES        :: EqualAI
        ES <--> EP      :: EqualAI
        EP ⇌ E + P        :: EqualAI
        E <--> F            :: OnlyA
        E + I ⇌ FI        :: OnlyA
        F + I ⇌ FI        :: EqualAI
        F + S ⇌ FS        :: EqualAI
        FS <--> FP      :: EqualAI
        FP ⇌ F + P        :: EqualAI
    end
end
_testhelper_tag_sets["BS2"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ ES        :: OnlyA
        ES <--> EP      :: EqualAI
        F + P ⇌ EP        :: EqualAI
        F ⇌ E               :: OnlyA
    end
end
_testhelper_tag_sets["BS3"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S <--> ES     :: OnlyA
        ES <--> EP      :: EqualAI
        F + P ⇌ EP        :: EqualAI
        F <--> E            :: OnlyA
    end
end
_testhelper_tag_sets["BS4"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S ⇌ ES            :: EqualAI
        ES <--> E + P         :: EqualAI
        E + I ⇌ EI            :: OnlyA
        EI + S <--> EIS   :: EqualAI
        EIS <--> FI + P   :: EqualAI
        FI ⇌ F + I            :: EqualAI
        F <--> E                :: OnlyA
    end
end
_testhelper_tag_sets["ADbridge"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S <--> ES                     :: EqualAI
        ES <--> E + P                     :: EqualAI
        (E + I ⇌ EI, ES + I ⇌ EIS)  :: OnlyA
        EI + S ⇌ EIS                  :: EqualAI
    end
end
_testhelper_tag_sets["ADbridge_noturnover"] =
    @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S <--> ES                     :: EqualAI
        ES <--> E + P                     :: OnlyA
        (E + I ⇌ EI, ES + I ⇌ EIS)  :: OnlyA
        EI + S ⇌ EIS                  :: EqualAI
    end
end
_testhelper_tag_sets["bridge_drain"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S ⇌ ES        :: EqualAI
        ES <--> F + P     :: EqualAI
        F <--> E            :: OnlyA
        E + I ⇌ EI        :: OnlyA
        F + I ⇌ FI        :: OnlyA
        FI ⇌ EI         :: EqualAI
    end
end
_testhelper_tag_sets["bridge_deadend"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    modifiers: I, J
    steps: begin
        E + S <--> ES         :: EqualAI
        ES <--> E + P         :: EqualAI
        E + J <--> EJ         :: EqualAI
        E + I ⇌ EI            :: OnlyA
        EJ + I ⇌ EIJ      :: OnlyA
        EI + J ⇌ EIJ      :: EqualAI
    end
end
_testhelper_tag_sets["R2funnel"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ ES        :: OnlyA
        ES <--> EP      :: EqualAI
        F + P ⇌ EP        :: EqualAI
        F <--> E            :: OnlyA
        F + S ⇌ FS        :: EqualAI
        FS <--> FP      :: EqualAI
        F + P ⇌ FP        :: EqualAI
        FP ⇌ EP         :: EqualAI
    end
end
_testhelper_tag_sets["R2funnel_SS"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S <--> ES     :: OnlyA
        ES <--> EP      :: EqualAI
        F + P ⇌ EP        :: EqualAI
        F <--> E            :: OnlyA
        F + S ⇌ FS        :: EqualAI
        FS <--> FP      :: EqualAI
        F + P ⇌ FP        :: EqualAI
        FP ⇌ EP         :: EqualAI
    end
end
_testhelper_tag_sets["S1"] = @_testhelper_unguarded @enzyme_mechanism begin
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
# Fuzz mechanism fz101_10232 (seed 101).
_testhelper_tag_sets["fz101_10232"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S <--> ES                             :: OnlyA
        ES <--> EP                              :: OnlyA
        EP <--> F + P                             :: EqualAI
        F <--> E                                    :: OnlyA
        (ES + I <--> EIS, F + I <--> FI)    :: OnlyA
        EP + S ⇌ ES + P                         :: EqualAI
    end
end

@testset "how an excluded step is switched off" begin
    D = EnzymeRates.Derivation
    rng = Random.MersenneTwister(2611)
    # The S binding takes up one metabolite into ES, which the inactive conformation's
    # P binding no longer joins to E: it loses its affinity. The chemistry takes up
    # nothing, so it has no affinity to lose; at steady state it is stopped.
    k_type = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES        :: OnlyA
            ES <--> EP      :: OnlyA
            E + P ⇌ EP        :: EqualAI
        end
    end
    off = D._built_reading(k_type)[:, 2]
    @test off[_testhelper_edge(k_type, "E + S ⇌ ES")] === D._LOSES_FORWARD
    @test off[_testhelper_edge(k_type, "ES <--> EP")] === D._STOPPED
    @test off[_testhelper_edge(k_type, "E + P ⇌ EP")] === D._RUNS
    @test all(==(D._RUNS), D._built_reading(k_type)[:, 1])
    @test D._guard(k_type).fault === D._NO_FAULT
    @test all(_testhelper_admission_agreement(k_type; rng))
    # With the chemistry untagged, ES stays, made from P: thermodynamics pins the S
    # binding's equilibrium constant, and a rapid-equilibrium step cannot turn slow.
    pinned = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES        :: OnlyA
            ES <--> EP      :: EqualAI
            E + P ⇌ EP        :: EqualAI
        end
    end
    v = D._guard(pinned)
    @test (v.fault, v.conformation, v.edge) ==
          (D._FAST_STEP_PINNED, 2, _testhelper_edge(pinned, "E + S ⇌ ES"))
    @test _testhelper_epath_verdict(pinned; rng) === :no_path
    # The same binding at steady state is stopped: S neither binds nor leaves in the
    # inactive conformation, and ES is made from P there.
    shut_door = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> ES     :: OnlyA
            ES <--> EP      :: EqualAI
            E + P ⇌ EP        :: EqualAI
        end
    end
    @test D._STOPPED in D._built_reading(shut_door)[:, 2]
    @test D._guard(shut_door).fault === D._NO_FAULT
    @test all(_testhelper_admission_agreement(shut_door; rng))
    # A rapid-equilibrium step that binds B and reacts at once, into EPQ, which the
    # P binding still makes: it cannot be switched off.
    fused_re = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA         :: EqualAI
            E + Q ⇌ EQ            :: EqualAI
            EA + B ⇌ EPQ      :: OnlyA
            EQ + P <--> EPQ   :: EqualAI
        end
    end
    @test D._guard(fused_re).fault === D._FAST_STEP_PINNED
    # At steady state the same step is stopped.
    fused_ss = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA         :: EqualAI
            E + Q ⇌ EQ            :: EqualAI
            EA + B <--> EPQ   :: OnlyA
            EQ + P <--> EPQ   :: EqualAI
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
    # steps the inactive conformation keeps do not join its two sides: EP lies in a group
    # of its own, and the isomerization loses its rate into EP, staying fast. The search
    # admits the mechanism, and the oracles agree on every face.
    D = EnzymeRates.Derivation
    am = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES        :: EqualAI
            ES ⇌ EP         :: OnlyA
            E + P <--> EP     :: OnlyA
        end
    end
    iso = _testhelper_edge(am, "ES ⇌ EP")
    @test D._admission_interior(am, D._built_reading(am)) ==
          D._Verdict(D._FAST_STEP_PINNED, D._LIMIT_UNTESTED, 0, 2, iso)
    @test D._guard(am) == D._WELL_POSED
    @test D._admitting_reading(am)[iso, 2] === D._LOSES_FORWARD
    @test all(_testhelper_admission_agreement(am; rng = Random.MersenneTwister(2615),
                                              corners = true))
    _testhelper_check_law_against_limit(am)
end

@testset "a group that loses one affinity names its step that cannot be switched off" begin
    # EA is a dead end, cut off in the inactive conformation, so the A binding at E
    # loses its affinity; its groupmate binds A into EPQ, which the P binding still
    # makes, at rapid equilibrium: the verdict names the groupmate.
    D = EnzymeRates.Derivation
    am = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EB + A ⇌ EPQ)  :: OnlyA
            E + B ⇌ EB                        :: EqualAI
            EQ + P <--> EPQ               :: EqualAI
            E + Q ⇌ EQ                        :: EqualAI
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
    # one metabolite from EA, cut off, into EAB, cut off too, so it loses its
    # affinity, and free E reaches both along the bindings.
    D = EnzymeRates.Derivation
    am = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA            :: OnlyA
            EA + B ⇌ EAB      :: OnlyA
            EAB <--> EPQ    :: OnlyA
            EPQ ⇌ EQ + P      :: OnlyA
            EQ ⇌ E + Q            :: OnlyA
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
    pinned = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES        :: OnlyA
            ES <--> EP      :: EqualAI
            E + P ⇌ EP        :: EqualAI
        end
        regulatory_site(2): R::OnlyI
    end
    shut_door = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S <--> ES     :: OnlyA
            ES <--> EP      :: EqualAI
            E + P ⇌ EP        :: EqualAI
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
    second_try_a = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            (E + D ⇌ ED, EBQ + D ⇌ EBDQ, EPQ + D ⇌ EDPQ,
             EQ + D ⇌ EDQ) :: EqualAI
            (E + Q ⇌ EQ, EBD + Q ⇌ EBDQ, ED + Q ⇌ EDQ,
             EDP + Q ⇌ EDPQ) :: OnlyA
            EBD + A <--> EPQR :: OnlyA
            (ED + B ⇌ EBD, EDQ + B ⇌ EBDQ, EQ + B ⇌ EBQ) :: OnlyA
            (ED + P <--> EDP, EDQ + P <--> EDPQ,
             EQ + P <--> EPQ, EQR + P <--> EPQR) :: EqualAI
            (EPQ + R ⇌ EPQR, EQ + R ⇌ EQR) :: EqualAI
        end
    end
    third_try = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            (E + D ⇌ ED, EAQ + D ⇌ EADQ, EPQ + D ⇌ EDPQ,
             EQ + D ⇌ EDQ) :: OnlyA
            (E + Q ⇌ EQ, EAD + Q ⇌ EADQ, ED + Q ⇌ EDQ,
             EDP + Q ⇌ EDPQ) :: EqualAI
            (EAD + B ⇌ EPQR, ED + B ⇌ EBD) :: OnlyA
            (EBD + A <--> EPQR, ED + A <--> EAD,
             EDQ + A <--> EADQ, EQ + A <--> EAQ) :: EqualAI
            (ED + P ⇌ EDP, EDQ + P ⇌ EDPQ, EQ + P ⇌ EPQ) :: EqualAI
            EPQ + R <--> EPQR :: OnlyA
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
        _testhelper_check_law_against_limit(am)
    end
    # The outward reading admits second_try_a but not third_try, so each later try is
    # needed.
    @test D._admission(second_try_a,
                       D._built_reading(second_try_a, D._OUTWARD_LOSE_AFFINITY)) ==
          D._WELL_POSED
    @test D._admission(third_try,
                       D._built_reading(third_try, D._OUTWARD_LOSE_AFFINITY)).fault ===
          D._BRIDGE
    # bridge_fast_a: its bridges, the rapid-equilibrium :OnlyA bindings such as the P
    # binding EDQ + P ⇌ EDPQ, lose their affinity under every reading, and stopping its
    # steady-state uptakes leaves forms that free E does not reach. The guard names the
    # first bridge in its search order, the B binding E_res_plusA_plusD_minusR + B ⇌ EPQ.
    # At zero P the ε-path oracle finds that the limit depends on the direction: the
    # rejection is real.
    bridge_fast_a =@_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            (E + D ⇌ ED, EAP + D ⇌ EADP, EP + D ⇌ EDP,
             EPQ + D ⇌ EDPQ) :: EqualAI
            (E + P ⇌ EP, ED + P ⇌ EDP, EDQ + P ⇌ EDPQ,
             ER_res_plusA_plusD_minusR + P ⇌ EADP) :: OnlyA
            E_res_plusA_plusD_minusR + B ⇌ EPQ :: OnlyA
            E_res_plusA_plusD_minusR + R ⇌ ER_res_plusA_plusD_minusR :: EqualAI
            (EAP + Q <--> EAPQ, ED + Q <--> EDQ,
             EDP + Q <--> EDPQ, EP + Q <--> EPQ,
             ER_res_plusA_plusD_minusR + Q <--> EADQ) :: OnlyA
            (ED + A <--> ER_res_plusA_plusD_minusR, EDP + A <--> EADP,
             EDQ + A <--> EADQ, EP + A <--> EAP,
             EPQ + A <--> EAPQ) :: EqualAI
        end
    end
    v = D._guard(bridge_fast_a)
    @test v == D._admission(bridge_fast_a, D._built_reading(bridge_fast_a))
    @test (v.fault, v.conformation, v.edge) ==
          (D._BRIDGE, 2,
           _testhelper_edge(bridge_fast_a, "E_res_plusA_plusD_minusR + B ⇌ EPQ"))
    @test D._admission(bridge_fast_a,
                       D._built_reading(bridge_fast_a, D._ALL_STOPPED)).fault === D._ISLAND
    @test D._admission(bridge_fast_a,
                       D._built_reading(bridge_fast_a, D._OUTWARD_LOSE_AFFINITY)).fault ===
          D._BRIDGE
    @test _testhelper_epath_verdict(bridge_fast_a; rng, face = Set([:P])) ===
          :direction_dependent
    # bridge_fast_b: every excluded uptake into a cut-off form is a rapid-equilibrium step,
    # so there is no steady-state one to read again.
    bridge_fast_b = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            (E + A ⇌ EA, EBQ + A ⇌ EABQ,
             EP_res_plusA_plusB_minusP + A ⇌ EAP_res_plusA_plusB_minusP,
             EQ + A ⇌ EAQ) :: OnlyA
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ,
             EP_res_plusA_plusB_minusP + Q ⇌ EABQ) :: OnlyA
            (E_res_plusA_plusB_minusP + P ⇌ EP_res_plusA_plusB_minusP,
             EA + P ⇌ EAP,
             EP_res_plusA_plusB_minusP + P ⇌ EABP) :: EqualAI
            E_res_plusA_plusB_minusP + D ⇌ EQR :: EqualAI
            (EA + B <--> EP_res_plusA_plusB_minusP, EAP + B <--> EABP,
             EAQ + B <--> EABQ,
             EP_res_plusA_plusB_minusP + B <--> EBP_res_plusA_plusB_minusP,
             EQ + B <--> EBQ) :: EqualAI
            (EQ + R <--> EQR,
             EQR + R <--> EDR_res_plusA_plusB_minusP) :: EqualAI
        end
    end
    @test D._guard(bridge_fast_b).fault === D._BRIDGE
    for c in 1:2
        @test D._built_reading(bridge_fast_b)[:, c] ==
              D._built_reading(bridge_fast_b, D._ALL_STOPPED)[:, c]
    end
    # The first try decides when it passes: with its S binding stopped, the K-type
    # uni-uni mechanism would leave ES out of free E's reach.
    k_type = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> ES     :: OnlyA
            ES <--> EP      :: OnlyA
            EP ⇌ E + P        :: EqualAI
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
    # both_readings, inactive conformation: the steps it keeps join EA to EAP and
    # EP to EPQ, each part one uptake from free E's. The B binding, from one part to
    # the other, must stop, or the steps that lose their affinity close a cycle with the
    # A binding at EP; the P binding at E must then lose its affinity, or no such step
    # leads to EP and EPQ.
    both_readings = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EP + A ⇌ EAP)  :: OnlyA
            E + P <--> EP                     :: OnlyA
            EA + P <--> EAP               :: EqualAI
            EA + B <--> EPQ               :: OnlyA
            EP + Q ⇌ EPQ                  :: EqualAI
        end
    end
    p = _testhelper_edge(both_readings, "E + P <--> EP")      # leads away: distance 0 to 1
    b = _testhelper_edge(both_readings, "EA + B <--> EPQ")    # distance 1 to 1
    a = _testhelper_edge(both_readings, "EP + A ⇌ EAP")       # rapid equilibrium, 1 to 1
    L, S = D._LOSES_FORWARD, D._STOPPED
    for (reading, name, at_p, at_b) in (
            (D._ALL_LOSE_AFFINITY, :all_lose_affinity, L, L),
            (D._OUTWARD_LOSE_AFFINITY, :outward_lose_affinity, L, S),
            (D._ALL_STOPPED, :all_stopped, S, S))
        off = D._built_reading(both_readings, reading)[:, 2]
        @test off[p] === at_p && off[b] === at_b
        @test off[a] === D._LOSES_FORWARD                     # a fast step is read one way
        @test _testhelper_reading(D._built_reading(both_readings, reading))[:, 2] ==
              _testhelper_reading(both_readings, name)[:, 2]
    end
    # The same need with the substrates' roles swapped and the P binding at EB :OnlyI.
    both_readings_only_i = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + B ⇌ EB, EP + B ⇌ EBP)  :: OnlyA
            E + P <--> EP                     :: OnlyA
            EB + P <--> EBP               :: OnlyI
            EB + A <--> EPQ               :: OnlyA
            EP + Q ⇌ EPQ                  :: EqualAI
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
        _testhelper_check_law_against_limit(am)
    end
    # outward_pumped: the first try admits it, and the outward reading alone would lose
    # it.
    outward_pumped = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            E + A ⇌ EA                                    :: OnlyA
            E + R <--> ER                                 :: OnlyA
            (EA + B <--> EAB, ER + B <--> EBR)  :: OnlyA
            (EA + P ⇌ EAP, ER + P ⇌ EPR)        :: EqualAI
            EAB + D <--> EPR + Q                    :: EqualAI
        end
    end
    @test D._admission(outward_pumped, D._built_reading(outward_pumped)) == D._WELL_POSED
    @test D._admission(outward_pumped,
                       D._built_reading(outward_pumped, D._OUTWARD_LOSE_AFFINITY)).fault ===
          D._PUMPED
    @test D._guard(outward_pumped) == D._WELL_POSED
    @test all(_testhelper_admission_agreement(outward_pumped;
                                              rng = Random.MersenneTwister(1)))
    # In the inactive conformation no uptake leads from free E to E_res_plusA_minusP,
    # which only the chemistry and the P release reach: its B binding leaves a part no
    # uptake reaches, so the outward reading stops it.
    unreached = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                            :: EqualAI
            EA <--> EP_res_plusA_minusP                        :: OnlyA
            EP_res_plusA_minusP ⇌ E_res_plusA_minusP + P      :: OnlyA
            E_res_plusA_minusP + B <--> EB_res_plusA_minusP   :: OnlyA
            EB_res_plusA_minusP <--> EQ                        :: OnlyA
            EQ ⇌ E + Q                                            :: OnlyA
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
    # A ping-pong with both chemistry steps :OnlyA. In the
    # inactive conformation the covalent forms F, FB and FP form one cut-off group
    # (F is E_res_plusA_minusP). EA <--> FP loses its stored direction, into that
    # group, and FB ⇌ EQ its reverse, also into that group: both point uphill from
    # free E's group, and the cut-off group empties through both halves. No built reading
    # loses a step in reverse.
    D = EnzymeRates.Derivation
    am = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                            :: EqualAI
            E + Q ⇌ EQ                                            :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP      :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP      :: EqualAI
            EA <--> EP_res_plusA_minusP                        :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                           :: OnlyA
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
    # FB ⇌ EQ in copy 2: the reverse (EQ → FB) vanishes, the forward does not,
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
    # A step keeps its kind, so a reading that stops the rapid-equilibrium FB ⇌ EQ is
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
    # guard's first built reading FB ⇌ EQ takes up nothing, so it cannot be switched
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
    # `filled`: a ping-pong whose B binding at the covalent form F (E_res_plusA_minusP)
    # and first chemistry step are :OnlyA. In the inactive conformation F and FP form a
    # cut-off group: EA <--> FP loses its stored direction, into it, and F + B ⇌ FB
    # its reverse, the release of B, also into it. With B zero the uptake of B dies, and
    # the release, fast though it vanishes, runs one way from the law's form FB into F,
    # whose group empties only through the slow chemistry FP → EA: F fills faster
    # than it empties.
    filled = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                            :: EqualAI
            E + Q ⇌ EQ                                            :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP      :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP      :: EqualAI
            EA <--> EP_res_plusA_minusP                        :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                           :: EqualAI
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
    # out, and offers the fix the guard admits: F + B <--> FB steady-state.
    @test D._steady_state_passes(filled, b_binding)
    @test D._admission_message(filled, v, r) ==
          "Mechanism: at zero A and B, in the inactive conformation the :OnlyA " *
          "rapid-equilibrium step E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP " *
          "runs one way into E_res_plusA_minusP, at a rate that vanishes yet grows " *
          "with the rapid-equilibrium speed, and the forms EP_res_plusA_minusP, " *
          "E_res_plusA_minusP leave only through the slower step " *
          "EA <--> EP_res_plusA_minusP, so E_res_plusA_minusP fills faster than " *
          "it empties and its weight depends on how fast the tagged step vanishes. Make " *
          "the kinetic group of E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP " *
          "steady-state."
    filled_ss = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                            :: EqualAI
            E + Q ⇌ EQ                                            :: EqualAI
            E_res_plusA_minusP + B <--> EB_res_plusA_minusP   :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP      :: EqualAI
            EA <--> EP_res_plusA_minusP                        :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                           :: EqualAI
        end
    end
    @test D._guard(filled_ss) == D._WELL_POSED
    @test all(_testhelper_admission_agreement(filled_ss; rng = Random.MersenneTwister(2653),
                                              corners = true))
    _testhelper_check_law_against_limit(filled_ss)
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
    pingpong = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                            :: EqualAI
            E + Q ⇌ EQ                                            :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP      :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP      :: EqualAI
            EA <--> EP_res_plusA_minusP                        :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                           :: OnlyA
        end
    end
    p = fill(D._RUNS, length(pingpong.catalytic.steps), 2)
    p[_testhelper_edge(pingpong, "EA <--> EP_res_plusA_minusP"), 2] = D._LOSES_FORWARD
    p[_testhelper_edge(pingpong, "EB_res_plusA_minusP ⇌ EQ"), 2] = D._LOSES_REVERSE
    @test D._admission(pingpong, p) == D._WELL_POSED

    # `slow_source`: a ping-pong whose B binding at F gives EQ at once, every step
    # :OnlyA but F + P ⇌ FP. In the inactive conformation free E stands alone and the
    # steps that leave it are steady-state uptakes, which lose their uptake. With B zero
    # F + B ⇌ EQ runs one way, fast though it vanishes, from EQ into F; but enzyme
    # reaches EQ only through steady-state steps, so the flow into F vanishes, and the
    # faces with B zero pass.
    slow_source = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A <--> EP_res_plusA_minusP, EQ + A <--> EAQ)    :: OnlyA
            (E + Q <--> EQ, EP_res_plusA_minusP + Q <--> EAQ)    :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP            :: EqualAI
            E_res_plusA_minusP + B ⇌ EQ                              :: OnlyA
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
    # At zero A and B, the release EQ → F + B is fast, vanishing and alive while the
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
    # F + B ⇌ EQ loses its affinity, and the direction it keeps is the release of B,
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
        v = D._face_condition(eps, law.present, faces[i], first(am.catalytic.substrates))
        v.fault === BR && D._admission_face(am, eps, law, faces, i) == D._WELL_POSED
    end
    # `set_aside`: a ping-pong whose chemistry steps and the bindings of A and of P are
    # :OnlyA. In the inactive conformation the law's forms are free E and EQ, and the
    # cut-off forms EA, FP, F and FB (F is E_res_plusA_minusP), which the excluded
    # steps join, are entered from free E by the A binding and from EQ by
    # FB ⇌ EQ, both losing their direction into them. Their vanishing steps open a
    # route from free E to EQ on which the reaction turns over, so the face condition
    # refuses the interior as a bridge. But the law's own rapid-equilibrium step
    # E + Q ⇌ EQ holds the two forms they touch in equilibrium, their fast classes drain,
    # and none fills faster than it empties: they are set aside.
    set_aside = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                            :: OnlyA
            E + Q ⇌ EQ                                            :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP      :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP      :: OnlyA
            EA <--> EP_res_plusA_minusP                        :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                           :: OnlyA
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
    # uptake E + A → EA is dead, so that pair runs one way, through its release, which
    # does not vanish and only empties EA: it is safe, and the cut-off forms go again.
    @test _testhelper_faces(set_aside)[4] == Set([:A])
    let eps = D._eps_network(set_aside, D._re_flags(set_aside.catalytic), r),
        law = D._law_network(set_aside), faces = D._faces(set_aside.catalytic)
        x = first(set_aside.catalytic.substrates)
        @test [D._face_condition(eps, law.present, f, x).fault for f in faces] ==
              [BR, NO, NO, BR, NO, NO, NO]
    end
    @test _testhelper_passes_only_set_aside(set_aside, r, 1)
    @test _testhelper_passes_only_set_aside(set_aside, r, 4)
    for face in _testhelper_faces(set_aside)
        @test _testhelper_epath_verdict(set_aside; rng, face, corners = true,
                                        reading = _testhelper_reading(r)) === :converges
    end
    # `kept`: the same ping-pong with E + Q <--> EQ steady-state, so free E and EQ are
    # no longer held in equilibrium, and EA ⇌ FP at rapid equilibrium, so the route
    # through the cut-off forms runs fast, its rates vanishing only as the reading says.
    # In the interior every transition at the cut-off forms runs both ways, and each of
    # them empties into the law's forms along fast steps that do not vanish, so none lies
    # in a fast class that could fail to drain; but the route joins two law forms whose
    # ratio the law sets otherwise: nothing is set aside, the interior stays a bridge, and
    # the limit there depends on how the two limits are taken.
    kept = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                            :: OnlyA
            E + Q <--> EQ                                         :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP      :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP      :: OnlyA
            EA ⇌ EP_res_plusA_minusP                           :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                           :: OnlyA
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
    # alone, and the cut-off forms EA, FP, F, FB and EQ touch nothing else.
    one_touch = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                            :: OnlyA
            E + Q ⇌ EQ                                            :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP      :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP      :: OnlyA
            EA <--> EP_res_plusA_minusP                        :: EqualAI
            EB_res_plusA_minusP ⇌ EQ                           :: EqualAI
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
    # The face condition alone refuses five faces. With P zero the uptake F + P → FP is
    # dead, and F + P ⇌ FP runs one way between two cut-off forms, through its release,
    # which vanishes and is fast; it is safe because no fast transition from the law's
    # forms reaches its source FP, so the flow through it vanishes.
    @test _testhelper_faces(one_touch)[6] == Set([:P])
    let eps = D._eps_network(one_touch, D._re_flags(one_touch.catalytic), o),
        law = D._law_network(one_touch), faces = D._faces(one_touch.catalytic)
        x = first(one_touch.catalytic.substrates)
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
    first_try = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            (E + A ⇌ EA, EDP + A ⇌ EADP, EP + A ⇌ EAP,
             EPR + A ⇌ EAPR) :: OnlyA
            (E + P ⇌ EP, EA + P ⇌ EAP, EAD + P ⇌ EADP,
             EAR + P ⇌ EAPR) :: OnlyA
            (EA + D ⇌ EAD, EAP + D ⇌ EADP, EAR + D ⇌ EADR,
             EP + D ⇌ EDP, EPR + D ⇌ EDPR) :: EqualAI
            (EA + R <--> EAR, EAD + R <--> EADR,
             EAP + R <--> EAPR, EDP + R <--> EDPR,
             EP + R <--> EPR) :: EqualAI
            EAD + B <--> EPQR :: OnlyA
            EPR + Q ⇌ EPQR :: OnlyA
        end
    end
    let reading = D._built_reading(first_try)
        eps = D._eps_network(first_try, D._re_flags(first_try.catalytic), reading)
        law, faces = D._law_network(first_try), D._faces(first_try.catalytic)
        x = first(first_try.catalytic.substrates)
        @test [D._face_condition(eps, law.present, f, x).fault for f in faces] ==
              [BR, NO, BR, BR, BR, BR, NO, NO, NO]
        @test D._admission(first_try, reading) == D._WELL_POSED
    end
    @test D._guard(first_try) == D._WELL_POSED
    @test all(_testhelper_admission_agreement(first_try; rng = Random.MersenneTwister(1)))
    _testhelper_check_law_against_limit(first_try)

    # R2funnel and R2funnel_SS under the guard's first built reading, at zero P. Their
    # cut-off forms touch free E alone, and each reaches it along transitions that do not
    # vanish; but F and FS form a fast class whose one slow way out, to FP, leads back
    # into it along fast steps, while the slow EP <--> ES leaves EP, which a fast
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
    # E + A <--> FP, F being E_res_plusA_minusP, rather than letting it lose its
    # uptake. In the inactive conformation free E stands alone, so the cut-off forms
    # touch one law form, and their fast classes drain through FP + Q <--> EAQ,
    # EAQ → EQ + A and EQ → E + Q. With A and B zero the stopped step runs one
    # way, from FP into free E, along its vanishing direction: a pair at a law form with
    # only its vanishing direction alive, which is not safe. It alone keeps them: nothing
    # is set aside, the face stays a bridge, and the limit there depends on the direction.
    stopped_exit = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A <--> EP_res_plusA_minusP, EQ + A <--> EAQ)    :: OnlyA
            (E + Q <--> EQ, EP_res_plusA_minusP + Q <--> EAQ)    :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP            :: EqualAI
            E_res_plusA_minusP + B ⇌ EQ                              :: OnlyA
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
        @test count(law.present[(length(stopped_exit.catalytic.enzyme_forms) + 1):end]) == 1
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
    # built reading. R-B1's cut-off forms EI and EIJ touch free E and the dead end
    # EJ, which only the steady-state E + J <--> EJ joins; R-B2's basins reach the
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
    # A six-step ping-pong with both chemistry steps :OnlyA (F is
    # E_res_plusA_minusP). No built reading switches off FB ⇌ EQ, a
    # rapid-equilibrium step that takes up nothing, so the guard searches. The first
    # reading it finds loses EA <--> FP's stored direction and FB ⇌ EQ's reverse,
    # both into the covalent group: the reading of "a reading may lose a step in either
    # direction".
    pingpong = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                            :: EqualAI
            E + Q ⇌ EQ                                            :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP      :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP      :: EqualAI
            EA <--> EP_res_plusA_minusP                        :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                           :: OnlyA
        end
    end
    pingpong_dimer = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                            :: EqualAI
            E + Q ⇌ EQ                                            :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP      :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP      :: EqualAI
            EA <--> EP_res_plusA_minusP                        :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                           :: OnlyA
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
        _testhelper_check_law_against_limit(am)
        @test all(_testhelper_admission_agreement(am; rng, corners = true))
    end
    # Its mirror image tags both chemistry steps :OnlyI, so the active conformation keeps
    # no chemistry: no reading helps.
    @test D._guard(_testhelper_mirror(pingpong)).fault === D._NO_ACTIVE_TURNOVER
    # The verdict does not depend on how the steps are written: the macro stores them in
    # canonical order. Nor on the order in which they are stored, which only orders the
    # search: stored backward, the same reading admits the mechanism.
    pingpong_rewritten = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            EB_res_plusA_minusP ⇌ EQ                           :: OnlyA
            EA <--> EP_res_plusA_minusP                        :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP      :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP      :: EqualAI
            E + Q ⇌ EQ                                            :: EqualAI
            E + A ⇌ EA                                            :: EqualAI
        end
    end
    @test pingpong_rewritten == pingpong
    @test D._guard(pingpong_rewritten) == D._WELL_POSED
    let cm = pingpong.catalytic, back = reverse(eachindex(cm.steps))
        backward = EnzymeRates.Mechanism(
            EnzymeRates.Mechanisms._unchecked_scheme(cm.metabolites, cm.enzyme_forms,
                                                     cm.substrates, cm.products,
                                                     cm.steps[back]),
            EnzymeRates.Mechanisms.Allostery(pingpong.allosteric.allosteric_states, 1,
                                             pingpong.allosteric.regulatory_sites))
        @test D._guard(backward) == D._WELL_POSED
        @test D._admitting_reading(backward) == r[back, :]
    end

    # The three-step ping-pong: the inactive conformation keeps only E + Q <--> EQ. The
    # built readings stop both chemistry steps, which leaves F out of free E's reach; the
    # first reading the search finds loses the first half's stored direction and the
    # second half's reverse, both into F.
    three_step = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> F + P        :: OnlyA
            E + Q <--> EQ         :: EqualAI
            F + B <--> EQ         :: OnlyA
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
    _testhelper_check_law_against_limit(three_step)
    @test all(_testhelper_admission_agreement(three_step; rng, corners = true))

    # A ping-pong with no limit under any reading: the covalent forms empty only by
    # binding B. The guard finds no reading, and the ε-path oracle converges on every face
    # under none of the consistent readings.
    no_limit = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EP_res_plusA_minusP                       :: OnlyA
            E + Q <--> EQ                                         :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP      :: EqualAI
            E_res_plusA_minusP + B ⇌ EQ                        :: OnlyA
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
    filled = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                            :: EqualAI
            E + Q ⇌ EQ                                            :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP      :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP      :: EqualAI
            EA <--> EP_res_plusA_minusP                        :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                           :: EqualAI
        end
    end
    @test D._admitting_reading(filled) === nothing
    @test D._guard(filled).fault !== D._NO_FAULT

    # A mechanism a built reading admits is admitted with it: `both_readings` of "a
    # mechanism may need one uptake to lose its affinity and another to stop".
    both_readings = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EP + A ⇌ EAP)  :: OnlyA
            E + P <--> EP                     :: OnlyA
            EA + P <--> EAP               :: EqualAI
            EA + B <--> EPQ               :: OnlyA
            EP + Q ⇌ EPQ                  :: EqualAI
        end
    end
    @test D._admitting_reading(both_readings) ==
          D._built_reading(both_readings, D._OUTWARD_LOSE_AFFINITY)

    # A rapid-equilibrium step inside one group is the fault no reading changes. The
    # message offers the two fixes, each of whose variants passes.
    pinned = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES        :: OnlyA
            ES <--> EP      :: EqualAI
            E + P ⇌ EP        :: EqualAI
        end
    end
    v = D._guard(pinned)
    @test v == D._Verdict(D._FAST_STEP_PINNED, D._LIMIT_UNTESTED, 0, 2,
                          _testhelper_edge(pinned, "E + S ⇌ ES"))
    @test D._guard_message(pinned) ==
          "Mechanism: in the inactive conformation the rapid-equilibrium step " *
          "E + S ⇌ ES is tagged :OnlyA, but the steps that conformation keeps still " *
          "join its two sides. Thermodynamics then fixes its equilibrium constant, and a " *
          "step keeps its kind in both conformations, so it cannot be switched off. Make " *
          "its kinetic group steady-state, or tag :OnlyA the steps that join its sides: " *
          "E + P ⇌ EP; ES <--> EP."
    @test _testhelper_consistent_readings(pinned) == Matrix{Symbol}[]

    # A cycle of rapid-equilibrium steps that turns the reaction over, each conformation
    # excluding two of them between its two groups. The built readings cannot switch off
    # the isomerizations, so under them the cycle is slow in each conformation; every
    # consistent reading lets each excluded step lose one direction, so the cycle runs
    # fast in both, and class 2 rejects it under the first, naming how each excluded step
    # of the cycle is switched off. The ε-path oracle finds the rate diverging.
    fast_cycle = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES        :: OnlyI
            ES ⇌ EstarS     :: OnlyA
            EstarS ⇌ EP     :: OnlyI
            E + P ⇌ EP        :: OnlyA
        end
    end
    @test D._verdict_under(fast_cycle, D._built_reading(fast_cycle)).fault ===
          D._FAST_STEP_PINNED
    v = D._guard(fast_cycle)
    @test v.fault === D._RE_TURNOVER
    reported = _testhelper_reading(D._guard_with_reading(fast_cycle).reading)
    @test _testhelper_reading_is_consistent(fast_cycle, reported)
    msg = D._guard_message(fast_cycle)
    @test occursin("(whose affinity vanishes in the ", msg)
    @test occursin("(whose rate one way vanishes in the ", msg)
    @test _testhelper_epath_verdict(fast_cycle; rng, reading = reported) === :diverges

    # A mechanism and its mirror image, every :OnlyA and :OnlyI swapped, get one verdict
    # when both conformations catalyze: the six-step ping-pong with a ternary-complex path
    # beside it, which both conformations keep. The search admits both, the mirror with
    # the same reading in the other conformation, and the oracles agree on each.
    pp_ternary = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                            :: EqualAI
            E + Q ⇌ EQ                                            :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP      :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP      :: EqualAI
            EA <--> EP_res_plusA_minusP                        :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                           :: OnlyA
            EA + B ⇌ EAB                                      :: EqualAI
            EAB <--> EPQ                                    :: EqualAI
            EQ + P ⇌ EPQ                                      :: EqualAI
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
        _testhelper_check_law_against_limit(am)
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
    large_admitted = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            E + D ⇌ ED :: EqualAI
            E + R ⇌ ER :: OnlyA
            (EAR + P <--> EAPR, EQR + P <--> EPQR,
             ER + P <--> EPR) :: OnlyA
            (EBD + A <--> EPQR, EPR + A <--> EAPR,
             ER + A <--> EAR) :: OnlyA
            (EBR + Q ⇌ EBQR, EPR + Q ⇌ EPQR,
             ER + Q ⇌ EQR) :: OnlyA
            (ED + B ⇌ EBD, EQR + B ⇌ EBQR, ER + B ⇌ EBR) :: OnlyA
        end
    end
    large_island = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            (E + B <--> EB, ED + B <--> EQ_res_plusB_plusD_minusQ) :: OnlyA
            (E + D <--> ED, EB + D <--> EQ_res_plusB_plusD_minusQ,
             EP + D <--> EDP,
             EQ_res_plusB_plusD_minusQ + D <--> EDQ_res_plusB_plusD_minusQ) :: OnlyA
            (E + P <--> EP, ED + P <--> EDP, ER + P <--> EPR) :: OnlyA
            (E + R <--> ER, EP + R <--> EPR) :: EqualAI
            E_res_plusB_plusD_minusQ + A ⇌ EPR :: OnlyA
            (E_res_plusB_plusD_minusQ + Q ⇌ EQ_res_plusB_plusD_minusQ,
             ED + Q ⇌ EDQ) :: OnlyA
        end
    end
    large_pumped = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            (E + A <--> EA, EB + A <--> EP_res_plusA_plusB_minusP,
             ER + A <--> EAR) :: OnlyA
            (E + B <--> EB, EA + B <--> EP_res_plusA_plusB_minusP) :: OnlyA
            (E + R <--> ER, EA + R <--> EAR,
             EQR + R <--> EDR_res_plusA_plusB_minusP) :: OnlyA
            (E_res_plusA_plusB_minusP + P ⇌ EP_res_plusA_plusB_minusP,
             EQR + P ⇌ EDP_res_plusA_plusB_minusP) :: EqualAI
            (E_res_plusA_plusB_minusP + D <--> EQR,
             EP_res_plusA_plusB_minusP + D <--> EDP_res_plusA_plusB_minusP,
             ER + D <--> EDR) :: OnlyA
            ER + Q ⇌ EQR :: OnlyA
        end
    end
    for am in (large_admitted, large_island, large_pumped)
        @test _testhelper_choice_count(am) > D._MAX_CHOICES_PER_CONFORMATION
    end
    @test D._admitting_reading(large_admitted) == D._built_reading(large_admitted)
    @test D._guard(large_island) ==
          D._Verdict(D._TOO_MANY_READINGS, D._LIMIT_UNTESTED, 0, 2, 0)
    @test D._guard_message(large_island) ==
          "Mechanism: the inactive conformation can switch off the tagged " *
          "steps in more than 20000 ways, too many for the guard to search, and none of " *
          "the three ways it tries first fits the rule for switching them off, so it " *
          "cannot show that the rate law is the limit."
    @test D._guard(large_pumped).fault === D._PUMPED
    @test D._guard(large_pumped) ==
          D._admission(large_pumped, D._built_reading(large_pumped))
    names = Dict(name => o for (o, name) in _testhelper_switch_off_names)
    LF, LR = :loses_forward, :loses_reverse
    for (am, column) in ((large_island, [LF, LF, LF, :runs, LF, LF, LF, LF, LF, :runs, LF,
                                         LF, LR, LR]),
                         (large_pumped, [LF, LF, LF, LF, LF, LF, LR, :runs, LF, LF, LF, LF,
                                         :runs, LR]))
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
    too_many_pairs = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            (E + A <--> EA, ED + A <--> EP_res_plusA_plusD_minusP,
             EQ + A <--> EAQ) :: OnlyI
            (E + D <--> ED, EA + D <--> EP_res_plusA_plusD_minusP,
             ER + D <--> EDR) :: OnlyI
            (E + Q <--> EQ, EA + Q <--> EAQ,
             EQR + Q <--> EBQ_res_plusA_plusD_minusP, ER + Q <--> EQR) :: OnlyI
            (E + R <--> ER, ED + R <--> EDR, EQ + R <--> EQR) :: OnlyA
            (E_res_plusA_plusD_minusP + P ⇌ EP_res_plusA_plusD_minusP,
             EQR + P ⇌ EBP_res_plusA_plusD_minusP) :: OnlyA
            (E_res_plusA_plusD_minusP + B ⇌ EQR,
             EP_res_plusA_plusD_minusP + B ⇌ EBP_res_plusA_plusD_minusP,
             EQ + B ⇌ EBQ) :: EqualAI
        end
    end
    too_many_active = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            (E + A ⇌ EA, EBR + A ⇌ EABR, ER + A ⇌ EAR) :: OnlyI
            (E + R <--> ER, EA + R <--> EAR,
             EQ_res_plusA_plusB_minusQ + R <--> EABR) :: OnlyI
            E_res_plusA_plusB_minusQ + D <--> EPR :: OnlyI
            E_res_plusA_plusB_minusQ + Q ⇌ EQ_res_plusA_plusB_minusQ :: OnlyI
            (EA + B <--> EQ_res_plusA_plusB_minusQ, EAR + B <--> EABR,
             ER + B <--> EBR) :: OnlyI
            ER + P ⇌ EPR :: OnlyI
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
        @test !EnzymeRates.Mechanisms.has_turnover_cycle(
            am.catalytic, EnzymeRates.Mechanisms.steps_in_active_conformation(am))
        @test D._guard(am) == D._Verdict(D._NO_ACTIVE_TURNOVER, D._LIMIT_UNTESTED, 0, 1, 1)
    end
    @test D._guard_message(too_many_pairs) ==
          "Mechanism: without its :OnlyI steps E + A <--> EA; " *
          "E + D <--> ED; E + Q <--> EQ; EA + Q <--> EAQ; " *
          "EA + D <--> EP_res_plusA_plusD_minusP; " *
          "ED + A <--> EP_res_plusA_plusD_minusP; EQ + A <--> EAQ; " *
          "EQR + Q <--> EBQ_res_plusA_plusD_minusP; ER + D <--> EDR; " *
          "ER + Q <--> EQR the active conformation turns no substrate into " *
          "product, but it must catalyze the reaction."
    @test D._guard_message(too_many_active) ==
          "Mechanism: without its :OnlyI steps E + A ⇌ EA; E + R <--> ER; " *
          "EA + R <--> EAR; EA + B <--> EQ_res_plusA_plusB_minusQ; " *
          "EAR + B <--> EABR; EBR + A ⇌ EABR; " *
          "EQ_res_plusA_plusB_minusQ + R <--> EABR; ER + A ⇌ EAR; " *
          "ER + B <--> EBR; ER + P ⇌ EPR; " *
          "E_res_plusA_plusB_minusQ + D <--> EPR; " *
          "E_res_plusA_plusB_minusQ + Q ⇌ EQ_res_plusA_plusB_minusQ the active " *
          "conformation turns no " *
          "substrate into product, but it must catalyze the reaction. Tag their kinetic " *
          "groups :EqualAI or :NonequalAI."

    # The same two ways with an active conformation that turns the reaction over, built
    # on random-order ter-ter: one full path of bindings, chemistry and releases kept by
    # both conformations; the other nine bindings, a second chemistry through EstarABD
    # and a second release order of Q and R tagged. In `cap_active` all three are :OnlyI,
    # which leaves the active conformation more than `_MAX_CHOICES_PER_CONFORMATION` ways;
    # in `cap_pairs` the bindings are :OnlyI and the rest :OnlyA, and the two
    # conformations' consistent columns pair up more than `_MAX_READING_PAIRS` times. Each
    # built reading leaves a group that no switched-off step reaches (EstarABD under
    # the two that lose affinity, EB or EQ when all stop), so none is consistent, and
    # the guard reports the cap.
    cap_active = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            E + A <--> EA                 :: EqualAI
            EA + B <--> EAB           :: EqualAI
            EAB + D <--> EABD     :: EqualAI
            E + B <--> EB                 :: OnlyI
            E + D <--> ED                 :: OnlyI
            EA + D <--> EAD           :: OnlyI
            EB + A <--> EAB           :: OnlyI
            EB + D <--> EBD           :: OnlyI
            ED + A <--> EAD           :: OnlyI
            ED + B <--> EBD           :: OnlyI
            EAD + B <--> EABD     :: OnlyI
            EBD + A <--> EABD     :: OnlyI
            EABD <--> EPQR      :: EqualAI
            EABD <--> EstarABD  :: OnlyI
            EstarABD <--> EPQR  :: OnlyI
            EPQR <--> EQR + P     :: EqualAI
            EQR <--> ER + Q           :: EqualAI
            ER <--> E + R                 :: EqualAI
            EQR <--> EQ + R           :: OnlyI
            EQ <--> E + Q                 :: OnlyI
        end
    end
    cap_pairs = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            E + A <--> EA                 :: EqualAI
            EA + B <--> EAB           :: EqualAI
            EAB + D <--> EABD     :: EqualAI
            E + B <--> EB                 :: OnlyI
            E + D <--> ED                 :: OnlyI
            EA + D <--> EAD           :: OnlyI
            EB + A <--> EAB           :: OnlyI
            EB + D <--> EBD           :: OnlyI
            ED + A <--> EAD           :: OnlyI
            ED + B <--> EBD           :: OnlyI
            EAD + B <--> EABD     :: OnlyI
            EBD + A <--> EABD     :: OnlyI
            EABD <--> EPQR      :: EqualAI
            EABD <--> EstarABD  :: OnlyA
            EstarABD <--> EPQR  :: OnlyA
            EPQR <--> EQR + P     :: EqualAI
            EQR <--> ER + Q           :: EqualAI
            ER <--> E + R                 :: EqualAI
            EQR <--> EQ + R           :: OnlyA
            EQ <--> E + Q                 :: OnlyA
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
        @test EnzymeRates.Mechanisms.has_turnover_cycle(
            am.catalytic, EnzymeRates.Mechanisms.steps_in_active_conformation(am))
        @test !any(r -> _testhelper_reading_is_consistent(
                       am, _testhelper_reading(am, r)), keys(_testhelper_guard_readings))
    end
    @test D._guard(cap_active) ==
          D._Verdict(D._TOO_MANY_READINGS, D._LIMIT_UNTESTED, 0, 1, 0)
    @test D._guard_message(cap_active) ==
          "Mechanism: the active conformation can switch off the tagged steps " *
          "in more than 20000 ways, too many for the guard to search, and none of the " *
          "three ways it tries first fits the rule for switching them off, so it cannot " *
          "show that the rate law is the limit."
    @test D._guard(cap_pairs) ==
          D._Verdict(D._TOO_MANY_READINGS, D._LIMIT_UNTESTED, 0, 0, 0)
    @test D._guard_message(cap_pairs) ==
          "Mechanism: the two conformations together can switch off the tagged " *
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
    # The uni-uni S binding :OnlyA: the chemistry and the P binding join E to ES in the
    # inactive conformation, so the S binding cannot lose its affinity there, and at rapid
    # equilibrium it cannot be stopped either.
    am = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES        :: OnlyA
            ES <--> EP      :: EqualAI
            E + P ⇌ EP        :: EqualAI
        end
    end
    m = am.catalytic
    bind, release, chem = (_testhelper_edge(am, t) for t in
                           ("E + S ⇌ ES", "E + P ⇌ EP", "ES <--> EP"))
    # The chemistry is stored substrate side first, so it runs toward products as stored.
    @test m.enzyme_forms[[m.steps[chem].from, m.steps[chem].to]] == [:ES, :EP]
    v = EnzymeRates.Derivation._guard(am)
    @test v == EnzymeRates.Derivation._Verdict(FP, U, 0, 2, bind)
    @test EnzymeRates.Derivation._admission_interior(
        am, EnzymeRates.Derivation._built_reading(am)) == v
    @test _testhelper_epath_verdict(am; rng) === :no_path
    line(k) = EnzymeRates.Mechanisms.step_text(m, k)
    @test EnzymeRates.Derivation._guard_message(am) ==
          _testhelper_pinned_cause(am, bind, 2) * " Make its kinetic group steady-state, " *
          "or tag :OnlyA the steps that join its sides: $(line(release)); $(line(chem))."
    # The fixes pass the guard: the steps that join its sides tagged, and the binding
    # steady-state ("how an excluded step is switched off"); so does the binding untagged,
    # which the message does not offer.
    untagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES        :: EqualAI
            ES <--> EP      :: EqualAI
            E + P ⇌ EP        :: EqualAI
        end
    end
    joins_tagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES        :: OnlyA
            ES <--> EP      :: OnlyA
            E + P ⇌ EP        :: OnlyA
        end
    end
    @test EnzymeRates.Derivation._guard(untagged) == EnzymeRates.Derivation._WELL_POSED
    @test EnzymeRates.Derivation._guard(joins_tagged) == EnzymeRates.Derivation._WELL_POSED
    # AD1 and AD11: a mixed :OnlyA group that binds B into EPQ and into EBQ, forms
    # the untagged steps still make. At steady state the group is stopped, and the law is
    # the limit on every face; at rapid equilibrium (AD1 RE) it cannot be switched off.
    # AD3c: a ::Mod copy's rapid-equilibrium binding whose sides untagged steps join.
    ad1_ss = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                    :: EqualAI
            (EA + B <--> EPQ, EQ + B <--> EBQ)  :: OnlyA
            EQ + P ⇌ EPQ                              :: EqualAI
            E + Q ⇌ EQ                                    :: EqualAI
            E + B ⇌ EB                                    :: EqualAI
            EB + Q ⇌ EBQ                              :: EqualAI
        end
    end
    ad1_re = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                :: EqualAI
            (EA + B ⇌ EPQ, EQ + B ⇌ EBQ)    :: OnlyA
            EQ + P <--> EPQ                       :: EqualAI
            E + Q ⇌ EQ                                :: EqualAI
            E + B ⇌ EB                                :: EqualAI
            EB + Q ⇌ EBQ                          :: EqualAI
        end
    end
    ad11 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                            :: EqualAI
            (EA + B <--> EPQ, EQ + B <--> EBQ,
             EP + B <--> EBP)                                 :: OnlyA
            EQ + P ⇌ EPQ                                      :: EqualAI
            E + Q ⇌ EQ                                            :: EqualAI
            E + P ⇌ EP                                            :: EqualAI
            E + B ⇌ EB                                            :: EqualAI
            EB + Q ⇌ EBQ                                      :: EqualAI
        end
    end
    ad3c = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES                    :: EqualAI
            ES <--> EP                  :: OnlyA
            E + P ⇌ EP                    :: EqualAI
            E + S ⇌ ESmod          :: OnlyA
            ES + S ⇌ ESSmod    :: EqualAI
            ESmod + S ⇌ ESSmod    :: EqualAI
        end
    end
    for x in (ad1_ss, ad11)
        @test EnzymeRates.Derivation._guard(x) == EnzymeRates.Derivation._WELL_POSED
        @test all(_testhelper_admission_agreement(x; rng))
    end
    for (x, text) in ((ad1_re, "EA + B ⇌ EPQ"), (ad3c, "E + S ⇌ ESmod"))
        @test EnzymeRates.Derivation._guard(x) ==
              EnzymeRates.Derivation._Verdict(FP, U, 0, 2, _testhelper_edge(x, text))
        @test _testhelper_exclusion(x, 2).failure === :fast
        @test _testhelper_epath_verdict(x; rng) === :no_path
    end
    # On AD3c the message's steps that join the sides run from the binding's unbound side
    # to its bound side, and that fix passes there. On AD1 RE both fixes pass.
    message(x) = EnzymeRates.Derivation._guard_message(x)
    joins = join((EnzymeRates.Mechanisms.step_text(ad3c.catalytic,
                                                   _testhelper_edge(ad3c, t))
                  for t in ("E + S ⇌ ES", "ES + S ⇌ ESSmod", "ESmod + S ⇌ ESSmod")),
                 "; ")
    @test endswith(message(ad3c), "Tag :OnlyA the steps that join its sides: $joins.")
    ad1_re_joins = join((EnzymeRates.Mechanisms.step_text(ad1_re.catalytic,
                                                         _testhelper_edge(ad1_re, t))
                         for t in ("E + A ⇌ EA", "E + Q ⇌ EQ", "EQ + P <--> EPQ")), "; ")
    @test endswith(message(ad1_re), "Make its kinetic group steady-state, or tag :OnlyA " *
                                    "the steps that join its sides: $ad1_re_joins.")
    ad1_untagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                    :: EqualAI
            (EA + B <--> EPQ, EQ + B <--> EBQ)  :: EqualAI
            EQ + P ⇌ EPQ                              :: EqualAI
            E + Q ⇌ EQ                                    :: EqualAI
            E + B ⇌ EB                                    :: EqualAI
            EB + Q ⇌ EBQ                              :: EqualAI
        end
    end
    ad1_joins_tagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                    :: EqualAI
            (EA + B <--> EPQ, EQ + B <--> EBQ)  :: OnlyA
            EQ + P ⇌ EPQ                              :: EqualAI
            E + Q ⇌ EQ                                    :: OnlyA
            E + B ⇌ EB                                    :: OnlyA
            EB + Q ⇌ EBQ                              :: OnlyA
        end
    end
    ad1_re_untagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                :: EqualAI
            (EA + B ⇌ EPQ, EQ + B ⇌ EBQ)    :: EqualAI
            EQ + P <--> EPQ                       :: EqualAI
            E + Q ⇌ EQ                                :: EqualAI
            E + B ⇌ EB                                :: EqualAI
            EB + Q ⇌ EBQ                          :: EqualAI
        end
    end
    ad1_re_joins_tagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                :: EqualAI
            (EA + B ⇌ EPQ, EQ + B ⇌ EBQ)    :: OnlyA
            EQ + P <--> EPQ                       :: EqualAI
            E + Q ⇌ EQ                                :: OnlyA
            E + B ⇌ EB                                :: OnlyA
            EB + Q ⇌ EBQ                          :: OnlyA
        end
    end
    ad3c_untagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES                    :: EqualAI
            ES <--> EP                  :: OnlyA
            E + P ⇌ EP                    :: EqualAI
            E + S ⇌ ESmod          :: EqualAI
            ES + S ⇌ ESSmod    :: EqualAI
            ESmod + S ⇌ ESSmod    :: EqualAI
        end
    end
    ad3c_joins_tagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES                    :: OnlyA
            ES <--> EP                  :: OnlyA
            E + P ⇌ EP                    :: EqualAI
            E + S ⇌ ESmod          :: OnlyA
            ES + S ⇌ ESSmod    :: OnlyA
            ESmod + S ⇌ ESSmod    :: OnlyA
        end
    end
    # Each variant of AD1, AD1 RE, AD3c and AD11 with the group untagged, or with the steps
    # that join the sides of its B or S binding tagged, passes the guard.
    ad11_untagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                            :: EqualAI
            (EA + B <--> EPQ, EQ + B <--> EBQ,
             EP + B <--> EBP)                                 :: EqualAI
            EQ + P ⇌ EPQ                                      :: EqualAI
            E + Q ⇌ EQ                                            :: EqualAI
            E + P ⇌ EP                                            :: EqualAI
            E + B ⇌ EB                                            :: EqualAI
            EB + Q ⇌ EBQ                                      :: EqualAI
        end
    end
    ad11_joins_tagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                            :: EqualAI
            (EA + B <--> EPQ, EQ + B <--> EBQ,
             EP + B <--> EBP)                                 :: OnlyA
            EQ + P ⇌ EPQ                                      :: EqualAI
            E + Q ⇌ EQ                                            :: OnlyA
            E + P ⇌ EP                                            :: EqualAI
            E + B ⇌ EB                                            :: OnlyA
            EB + Q ⇌ EBQ                                      :: OnlyA
        end
    end
    for y in (ad1_untagged, ad1_joins_tagged, ad1_re_untagged, ad1_re_joins_tagged,
              ad3c_untagged, ad3c_joins_tagged, ad11_untagged, ad11_joins_tagged)
        @test EnzymeRates.Derivation._guard(y) == EnzymeRates.Derivation._WELL_POSED
    end
    # The message prints a fix only when its variant passes. A dead-end I joins the S
    # binding's sides a second way, through EI and EIS, so tagging the steps of the
    # path the message would name leaves that binding inside one component, and the guard
    # rejects its steady-state variant too: the message offers no fix. (Untagging the
    # binding passes, but is not one of the two fixes the message offers.)
    deadend = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES        :: OnlyA
            ES <--> EP      :: EqualAI
            E + P ⇌ EP        :: EqualAI
            E + I ⇌ EI        :: EqualAI
            ES + I ⇌ EIS  :: EqualAI
            EI + S ⇌ EIS  :: EqualAI
        end
    end
    deadend_untagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES        :: EqualAI
            ES <--> EP      :: EqualAI
            E + P ⇌ EP        :: EqualAI
            E + I ⇌ EI        :: EqualAI
            ES + I ⇌ EIS  :: EqualAI
            EI + S ⇌ EIS  :: EqualAI
        end
    end
    deadend_joins_tagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES        :: OnlyA
            ES <--> EP      :: OnlyA
            E + P ⇌ EP        :: OnlyA
            E + I ⇌ EI        :: EqualAI
            ES + I ⇌ EIS  :: EqualAI
            EI + S ⇌ EIS  :: EqualAI
        end
    end
    v = EnzymeRates.Derivation._guard(deadend)
    @test v == EnzymeRates.Derivation._Verdict(FP, U, 0, 2,
                                               _testhelper_edge(deadend, "E + S ⇌ ES"))
    @test _testhelper_epath_verdict(deadend; rng) === :no_path
    @test EnzymeRates.Derivation._guard_message(deadend) ==
          _testhelper_pinned_cause(deadend, v.edge, 2)
    @test EnzymeRates.Derivation._guard(deadend_untagged) ==
          EnzymeRates.Derivation._WELL_POSED
    @test EnzymeRates.Derivation._guard(deadend_joins_tagged).fault === FP
    # I's binding at E :OnlyA too: untagging it leaves the S binding inside one component,
    # while tagging the steps that join its sides separates every binding, so the message
    # offers only the joins fix.
    two_bindings = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES        :: OnlyA
            ES <--> EP      :: EqualAI
            E + P ⇌ EP        :: EqualAI
            E + I ⇌ EI        :: OnlyA
            ES + I ⇌ EIS  :: EqualAI
            EI + S ⇌ EIS  :: EqualAI
        end
    end
    two_bindings_untagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES        :: OnlyA
            ES <--> EP      :: EqualAI
            E + P ⇌ EP        :: EqualAI
            E + I ⇌ EI        :: EqualAI
            ES + I ⇌ EIS  :: EqualAI
            EI + S ⇌ EIS  :: EqualAI
        end
    end
    two_bindings_joins_tagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES        :: OnlyA
            ES <--> EP      :: OnlyA
            E + P ⇌ EP        :: OnlyA
            E + I ⇌ EI        :: OnlyA
            ES + I ⇌ EIS  :: OnlyA
            EI + S ⇌ EIS  :: OnlyA
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
    @test EnzymeRates.Derivation._guard_message(two_bindings) ==
          _testhelper_pinned_cause(two_bindings, v.edge, 2) *
          " Tag :OnlyA the steps that join its sides: $joins."
    @test EnzymeRates.Derivation._guard(two_bindings_untagged).fault === FP
    @test EnzymeRates.Derivation._guard(two_bindings_joins_tagged) ==
          EnzymeRates.Derivation._WELL_POSED
end

@testset "a tagged binding pinned inside one group, or closing a cycle of a reading" begin
    rng = Random.MersenneTwister(37)
    FP, U = EnzymeRates.Derivation._FAST_STEP_PINNED, EnzymeRates.Derivation._LIMIT_UNTESTED
    # The ter cube of seven :OnlyA edges: the untagged EA + C, EC + B and EB + A
    # bindings join EA to EAB, so its rapid-equilibrium :OnlyA binding of B stays
    # inside one component and cannot be switched off.
    cube = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, C
        products: P
        steps: begin
            E + A ⇌ EA             :: OnlyA
            E + B ⇌ EB             :: OnlyA
            E + C ⇌ EC             :: OnlyA
            EA + B ⇌ EAB       :: OnlyA
            EA + C ⇌ EAC       :: EqualAI
            EB + A ⇌ EAB       :: EqualAI
            EB + C ⇌ EBC       :: EqualAI
            EC + A ⇌ EAC       :: EqualAI
            EC + B ⇌ EBC       :: EqualAI
            EAB + C ⇌ EABC :: OnlyA
            EAC + B ⇌ EABC :: OnlyA
            EBC + A ⇌ EABC :: OnlyA
            EABC <--> EP     :: OnlyA
            E + P ⇌ EP             :: EqualAI
        end
    end
    @test _testhelper_exclusion(cube, 2).failure === :fast
    @test EnzymeRates.Derivation._guard(cube) ==
          EnzymeRates.Derivation._Verdict(FP, U, 0, 2,
                                          _testhelper_edge(cube, "EA + B ⇌ EAB"))
    @test _testhelper_epath_verdict(cube; rng) === :no_path
    # The message offers three fixes, each passing the guard: the B binding at EA
    # untagged, the five untagged bindings that join its sides (EA + C, EC + A,
    # EC + B, EB + C, EB + A) tagged, and the B binding steady-state.
    cube_untagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, C
        products: P
        steps: begin
            E + A ⇌ EA             :: OnlyA
            E + B ⇌ EB             :: OnlyA
            E + C ⇌ EC             :: OnlyA
            EA + B ⇌ EAB       :: EqualAI
            EA + C ⇌ EAC       :: EqualAI
            EB + A ⇌ EAB       :: EqualAI
            EB + C ⇌ EBC       :: EqualAI
            EC + A ⇌ EAC       :: EqualAI
            EC + B ⇌ EBC       :: EqualAI
            EAB + C ⇌ EABC :: OnlyA
            EAC + B ⇌ EABC :: OnlyA
            EBC + A ⇌ EABC :: OnlyA
            EABC <--> EP     :: OnlyA
            E + P ⇌ EP             :: EqualAI
        end
    end
    cube_joins_tagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, C
        products: P
        steps: begin
            E + A ⇌ EA             :: OnlyA
            E + B ⇌ EB             :: OnlyA
            E + C ⇌ EC             :: OnlyA
            EA + B ⇌ EAB       :: OnlyA
            EA + C ⇌ EAC       :: OnlyA
            EB + A ⇌ EAB       :: OnlyA
            EB + C ⇌ EBC       :: OnlyA
            EC + A ⇌ EAC       :: OnlyA
            EC + B ⇌ EBC       :: OnlyA
            EAB + C ⇌ EABC :: OnlyA
            EAC + B ⇌ EABC :: OnlyA
            EBC + A ⇌ EABC :: OnlyA
            EABC <--> EP     :: OnlyA
            E + P ⇌ EP             :: EqualAI
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
        EnzymeRates.Derivation._guard_message(cube),
        "Make its kinetic group steady-state, or tag :OnlyA the steps that join its " *
        "sides: $joins.")
    # Both bindings of the ordered bi-uni :OnlyA: the chemistry and the P binding join
    # EAB back to E, so EA alone is cut off. No built reading switches off the B
    # binding, which takes up B into a form free E's group holds; the one consistent
    # reading loses the A binding's uptake and the B binding's release, both into EA.
    # Their vanishing directions then bridge E and EAB, and the limit depends on the
    # direction.
    am = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            E + A ⇌ EA            :: OnlyA
            EA + B ⇌ EAB      :: OnlyA
            EAB <--> EP       :: EqualAI
            E + P ⇌ EP            :: EqualAI
        end
    end
    b = _testhelper_edge(am, "EA + B ⇌ EAB")
    @test _testhelper_exclusion(am, 2).failure === :fast
    @test EnzymeRates.Derivation.exclusion_components(am, :inactive) ==
          [1, 2, 1, 1]  # E, EA, EAB, EP
    @test EnzymeRates.Derivation.exclusion_components(am, :active) == [1, 1, 1, 1]
    v = EnzymeRates.Derivation._guard(am)
    @test v == EnzymeRates.Derivation._Verdict(EnzymeRates.Derivation._BRIDGE, U, 1, 2, b)
    @test _testhelper_epath_verdict(am; rng) === :no_path
    reported = _testhelper_reading(EnzymeRates.Derivation._guard_with_reading(am).reading)
    @test _testhelper_reading_is_consistent(am, reported)
    @test _testhelper_epath_verdict(am; rng, corners = true, reading = reported) ===
          :direction_dependent
    @test EnzymeRates.Derivation._guard_message(am) ==
          "Mechanism: the :OnlyA step EA + B ⇌ EAB bridges forms of the " *
          "inactive conformation that its other steps keep apart, so admission cannot " *
          "show that the limit does not depend on how fast that step vanishes. Make the " *
          "kinetic group of EA + B ⇌ EAB steady-state."
    # Both bindings untagged, or the chemistry and the P binding tagged too, the mechanism
    # passes the guard; untagging the B binding alone leaves the A binding's sides joined.
    untagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            E + A ⇌ EA            :: EqualAI
            EA + B ⇌ EAB      :: EqualAI
            EAB <--> EP       :: EqualAI
            E + P ⇌ EP            :: EqualAI
        end
    end
    joins_tagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            E + A ⇌ EA            :: OnlyA
            EA + B ⇌ EAB      :: OnlyA
            EAB <--> EP       :: OnlyA
            E + P ⇌ EP            :: OnlyA
        end
    end
    one_untagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            E + A ⇌ EA            :: OnlyA
            EA + B ⇌ EAB      :: EqualAI
            EAB <--> EP       :: EqualAI
            E + P ⇌ EP            :: EqualAI
        end
    end
    @test EnzymeRates.Derivation._guard(untagged) == EnzymeRates.Derivation._WELL_POSED
    @test EnzymeRates.Derivation._guard(joins_tagged) == EnzymeRates.Derivation._WELL_POSED
    @test EnzymeRates.Derivation._guard(one_untagged).fault === FP
    # With I forms, EA is still cut off alone, and the one consistent reading loses both
    # bindings into it, as above. Admission passes the faces with I present, and class 3
    # rejects zero I, where the law splits a rapid-equilibrium segment, as the RE-speed
    # oracle finds. Tagging the chemistry and the P binding too leaves the law a second
    # closed class on a face.
    with_i = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P
        modifiers: I
        steps: begin
            E + A ⇌ EA                :: OnlyA
            EA + B ⇌ EAB          :: OnlyA
            EAB <--> EP           :: EqualAI
            E + P ⇌ EP                :: EqualAI
            E + I ⇌ EI                :: EqualAI
            EAB + I ⇌ EABI    :: EqualAI
            EI + A ⇌ EAI          :: EqualAI
            EAI + B ⇌ EABI    :: EqualAI
        end
    end
    with_i_untagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P
        modifiers: I
        steps: begin
            E + A ⇌ EA                :: EqualAI
            EA + B ⇌ EAB          :: EqualAI
            EAB <--> EP           :: EqualAI
            E + P ⇌ EP                :: EqualAI
            E + I ⇌ EI                :: EqualAI
            EAB + I ⇌ EABI    :: EqualAI
            EI + A ⇌ EAI          :: EqualAI
            EAI + B ⇌ EABI    :: EqualAI
        end
    end
    with_i_joins_tagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P
        modifiers: I
        steps: begin
            E + A ⇌ EA                :: OnlyA
            EA + B ⇌ EAB          :: OnlyA
            EAB <--> EP           :: OnlyA
            E + P ⇌ EP                :: OnlyA
            E + I ⇌ EI                :: EqualAI
            EAB + I ⇌ EABI    :: EqualAI
            EI + A ⇌ EAI          :: EqualAI
            EAI + B ⇌ EABI    :: EqualAI
        end
    end
    v = EnzymeRates.Derivation._guard(with_i)
    @test (v.fault, _testhelper_faces(with_i)[v.face]) ==
          (EnzymeRates.Derivation._SPLIT, Set([:I]))
    @test _testhelper_epath_verdict(with_i; rng) === :no_path
    @test _testhelper_agrees(v.fault, _testhelper_respeed_verdicts(
        with_i; rng = Random.MersenneTwister(38), law = :limit)[v.face])
    @test EnzymeRates.Derivation._guard_message(with_i) ==
          "Mechanism: at zero I, the forms {E (active), EA (active), " *
          "EAB (active), EP (active), E (inactive), EP (inactive)} and " *
          "{EAB (inactive)} of one rapid-equilibrium segment are joined only by " *
          "steps that take up a zeroed metabolite (E + I ⇌ EI and " *
          "EAB + I ⇌ EABI), so the rate law cannot weight one against the other."
    @test EnzymeRates.Derivation._guard(with_i_untagged) ==
          EnzymeRates.Derivation._WELL_POSED
    @test EnzymeRates.Derivation._guard(with_i_joins_tagged).fault ===
          EnzymeRates.Derivation._NONUNIQUE
    # Under the guard's built readings the arcs close a cycle between two cut-off parts:
    # EA with EAQ, and EQ with EAB and EPQ. The B binding at EA enters
    # the second, the A binding at EQ the first, so the two affinities cannot vanish
    # together. A reading in which one of them loses its release instead has no cycle,
    # and the search admits the mechanism under it, the oracles agreeing.
    cycle = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA              :: OnlyA
            EA + Q ⇌ EAQ        :: EqualAI
            EA + B ⇌ EAB        :: OnlyA
            EAB <--> EPQ      :: EqualAI
            EQ + P ⇌ EPQ        :: EqualAI
            EQ + A ⇌ EAQ        :: OnlyA
            E + Q ⇌ EQ              :: OnlyA
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
    _testhelper_check_law_against_limit(cycle)
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
        _testhelper_check_law_against_limit(x)
    end
    _testhelper_edge_between(x, a, b) = only(k for (k, e) in enumerate(x.catalytic.steps)
        if Set(String.(x.catalytic.enzyme_forms[[e.from, e.to]])) == Set((a, b)))
    _testhelper_island_verdict(x) =
        _testhelper_epath_verdict(x; rng) in (:direction_dependent, :wrong_limit)
    # Both half-reactions of the ping-pong :OnlyA. Under the first built reading half 2
    # takes up Q into EB_res_plusA_minusP, which the inactive conformation cuts off, so
    # it loses its affinity and free E reaches the covalent forms along it: no island. At
    # zero A and B, B's binding feeds E_res_plusA_minusP, whose only way out, half 1,
    # is stopped: a pumped basin, as the ε-path oracle finds under that reading. The
    # search admits the ping-pong under a reading in which half 1 loses its rate into the
    # covalent forms instead.
    both = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA                                       :: EqualAI
            EA <--> E_res_plusA_minusP + P                   :: OnlyA
            E_res_plusA_minusP + B <--> EB_res_plusA_minusP :: EqualAI
            EB_res_plusA_minusP <--> E + Q                     :: OnlyA
        end
    end
    half1 = _testhelper_edge(both, "E_res_plusA_minusP + P <--> EA")
    feed = _testhelper_edge(both, "E_res_plusA_minusP + B <--> EB_res_plusA_minusP")
    @test both.catalytic.enzyme_forms[[both.catalytic.steps[half1].from,
                                       both.catalytic.steps[half1].to]] ==
          [:E_res_plusA_minusP, :EA]
    built = D._built_reading(both)
    v = D._admission(both, built)
    @test v == EnzymeRates.Derivation._Verdict(PU, U, 3, 2, feed)
    @test isempty(_testhelper_exclusion(both, 2).lost)
    @test _testhelper_admission_agrees(
        PU, _testhelper_epath_verdict(both; rng, face = _testhelper_faces(both)[3]))
    @test D._admission_message(both, v, built) ==
          "Mechanism: at zero A and B, " *
          "$(EnzymeRates.Mechanisms.step_text(both.catalytic, feed)) feeds the forms " *
          "E_res_plusA_minusP of the inactive conformation, which leave only through " *
          "the vanishing step $(EnzymeRates.Mechanisms.step_text(both.catalytic, half1)) " *
          "while nothing runs back the way it came, so admission cannot show that their " *
          "weight does not depend on how fast that step vanishes."
    _testhelper_admitted_by_search(both)
    @test EnzymeRates.Derivation._guard(_testhelper_tag_sets["pp_V_half2"]) ==
          EnzymeRates.Derivation._WELL_POSED
    # R3 and R4: a second free form behind a tagged isomerization; untagging it passes.
    # Under the built readings R3's isomerization is steady-state, so it is stopped and
    # F, FI are an island; R4's is at rapid equilibrium, which they cannot switch off.
    # Each isomerization joins F's group to free E's, and the search lets it lose its rate
    # into F.
    r3 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES        :: EqualAI
            ES <--> EP      :: EqualAI
            E + P ⇌ EP        :: EqualAI
            E <--> F            :: OnlyA
            F + I ⇌ FI        :: EqualAI
        end
    end
    r3_untagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES        :: EqualAI
            ES <--> EP      :: EqualAI
            E + P ⇌ EP        :: EqualAI
            E <--> F            :: EqualAI
            F + I ⇌ FI        :: EqualAI
        end
    end
    r4 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES        :: EqualAI
            ES <--> EP      :: NonequalAI
            E + P ⇌ EP        :: EqualAI
            F ⇌ E               :: OnlyA
        end
    end
    r4_untagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES        :: EqualAI
            ES <--> EP      :: NonequalAI
            E + P ⇌ EP        :: EqualAI
            F ⇌ E               :: EqualAI
        end
    end
    iso = _testhelper_edge_between(r3, "E", "F")
    # An isomerization of free forms ties on both sides; it is stored from E.
    @test r3.catalytic.steps[iso].from == 1
    @test D._admission_interior(r3, D._built_reading(r3)) ==
          EnzymeRates.Derivation._Verdict(IS, U, 0, 2, iso)
    @test _testhelper_island_verdict(r3)
    _testhelper_admitted_by_search(r3)
    iso = _testhelper_edge_between(r4, "E", "F")
    @test r4.catalytic.steps[iso].from == 1
    @test D._admission_interior(r4, D._built_reading(r4)) ==
          EnzymeRates.Derivation._Verdict(FP, U, 0, 2, iso)
    @test _testhelper_epath_verdict(r4; rng) === :no_path
    _testhelper_admitted_by_search(r4)
    for fixed in (r3_untagged, r4_untagged)
        @test EnzymeRates.Derivation._guard(fixed) == EnzymeRates.Derivation._WELL_POSED
    end
    # AD2b's P binding at F is at rapid equilibrium and joins F to EP, in free E's group,
    # so no built reading switches it off. Its consistent readings let it lose its release
    # into F; with P zero its uptake then dies, and it runs one way into F, which empties
    # only through the slow F <--> E: F fills faster than it empties, as the ε-path oracle
    # finds, and the message offers the steady-state fix, which the search admits. In the
    # ping-pong with B's binding and half 1 :OnlyA, the covalent form is an island behind
    # half 1 under the built readings; untagged, half 1 joins it to E, and B's steady-state
    # binding, its sides joined, is stopped, a fix that passes every face, the oracles
    # agreeing.
    ad2b = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES        :: EqualAI
            ES <--> EP      :: EqualAI
            F + P ⇌ EP        :: OnlyA
            F <--> E            :: OnlyA
        end
    end
    bbind = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA                                       :: EqualAI
            EA <--> E_res_plusA_minusP + P                   :: OnlyA
            E_res_plusA_minusP + B <--> EB_res_plusA_minusP :: OnlyA
            EB_res_plusA_minusP <--> E + Q                     :: EqualAI
        end
    end
    bbind_untagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA                                       :: EqualAI
            EA <--> E_res_plusA_minusP + P                   :: EqualAI
            E_res_plusA_minusP + B <--> EB_res_plusA_minusP :: OnlyA
            EB_res_plusA_minusP <--> E + Q                     :: EqualAI
        end
    end
    ad2b_ss = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES        :: EqualAI
            ES <--> EP      :: EqualAI
            F + P <--> EP     :: OnlyA
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
    @test EnzymeRates.Derivation._guard_message(ad2b) ==
          "Mechanism: at zero P, in the inactive conformation the :OnlyA " *
          "rapid-equilibrium step F + P ⇌ EP runs one way into F, at a rate that " *
          "vanishes yet grows with the rapid-equilibrium speed, and F leaves only " *
          "through the slower step F <--> E, so F fills faster than it empties and its " *
          "weight depends on how fast the tagged step vanishes. Make the kinetic group " *
          "of F + P ⇌ EP steady-state."
    _testhelper_admitted_by_search(ad2b_ss)
    half1 = _testhelper_edge(bbind, "E_res_plusA_minusP + P <--> EA")
    @test D._admission_interior(bbind, D._built_reading(bbind)) ==
          EnzymeRates.Derivation._Verdict(IS, U, 0, 2, half1)
    @test _testhelper_island_verdict(bbind)
    _testhelper_admitted_by_search(bbind)
    @test EnzymeRates.Derivation._guard(bbind_untagged) ==
          EnzymeRates.Derivation._WELL_POSED
    @test all(_testhelper_admission_agreement(bbind_untagged; rng))
    # Half 1 and its copy at EAI share one :OnlyA group, and I's binding at E is :OnlyA
    # too. As in the ping-pong above, under the first built reading, at zero A and B, the
    # covalent forms are fed, by B's binding and by the copy of half 1 at EAI running
    # forward, and their only ways out, the two copies of half 1 taking P back up, vanish.
    # The guard names the feed whose reverse is dead, B's binding: nothing runs back
    # through it. The search admits it as it admits that ping-pong. Untagging that group
    # joins EI to E, so I's rapid-equilibrium binding cannot be switched off.
    shared_half1 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            E + A <--> EA                                       :: EqualAI
            (EA <--> E_res_plusA_minusP + P,
             EAI <--> EI_res_plusA_minusP + P)             :: OnlyA
            E_res_plusA_minusP + B <--> EB_res_plusA_minusP :: EqualAI
            EB_res_plusA_minusP <--> E + Q                     :: OnlyA
            E + I ⇌ EI                                          :: OnlyA
            EI + A <--> EAI                                 :: EqualAI
            E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP    :: EqualAI
        end
    end
    shared_half1_untagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            E + A <--> EA                                       :: EqualAI
            (EA <--> E_res_plusA_minusP + P,
             EAI <--> EI_res_plusA_minusP + P)             :: EqualAI
            E_res_plusA_minusP + B <--> EB_res_plusA_minusP :: EqualAI
            EB_res_plusA_minusP <--> E + Q                     :: OnlyA
            E + I ⇌ EI                                          :: OnlyA
            EI + A <--> EAI                                 :: EqualAI
            E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP    :: EqualAI
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
          "Mechanism: at zero A and B, " *
          "$(shared_line("E_res_plusA_minusP + B <--> EB_res_plusA_minusP")) feeds the " *
          "forms EI_res_plusA_minusP, E_res_plusA_minusP of the inactive " *
          "conformation, which leave only through the vanishing steps " *
          "$(shared_line("EI_res_plusA_minusP + P <--> EAI")); " *
          "$(shared_line("E_res_plusA_minusP + P <--> EA")) while nothing runs back " *
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
    # B's fused binding into EPQ and its plain binding at EQ share one :OnlyA group;
    # in the inactive conformation no untagged step joins EQ to EBQ, so the plain
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
    # vanish (F <--> E, stopped) or take up P (F + P ⇌ EP), and EP feeds it by a
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
          "Mechanism: at zero P, $(line(k)) feeds the forms F of the inactive " *
          "conformation, which leave only through the vanishing step $(line(iso)) while " *
          "nothing runs back the way it came, so admission cannot show that their weight " *
          "does not depend on how fast that step vanishes."
    @test EnzymeRates.Derivation._guard(am) == EnzymeRates.Derivation._WELL_POSED
    @test EnzymeRates.Derivation._admitting_reading(am)[iso, 2] !==
          EnzymeRates.Derivation._STOPPED
    @test all(_testhelper_admission_agreement(am; rng, corners = true))
    _testhelper_check_law_against_limit(am)
end

@testset "a stranded basin with two vanishing exits names both steps in its message" begin
    # Under the guard's first built reading, at P zero, ES feeds F by a release whose
    # reverse takes up P; F and FI leave the inactive conformation's basin only through
    # the :OnlyA steps EI ⇌ FI and F <--> E, which vanish there, so admission's face
    # condition rejects that face. But that reading cannot switch off EI ⇌ FI, at
    # rapid equilibrium, which consistent readings let lose its rate into FI; under the
    # first of them class 3 rejects zero I, where the law splits a rapid-equilibrium
    # segment, as the RE-speed oracle finds.
    am = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + I ⇌ EI        :: NonequalAI
            E + S ⇌ ES        :: OnlyA
            EI ⇌ FI         :: OnlyA
            F <--> E            :: OnlyA
            F + P <--> ES     :: NonequalAI
            F + I ⇌ FI        :: OnlyA
        end
    end
    rng = Random.MersenneTwister(43)
    k = _testhelper_edge(am, "FI ⇌ EI")
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
    @test EnzymeRates.Derivation._guard_message(am) ==
          "Mechanism: at zero I, the forms {E (active), ES (active), " *
          "E (inactive)} and {F (active)} of one rapid-equilibrium segment are joined " *
          "only by steps that take up a zeroed metabolite (E + I ⇌ EI and " *
          "F + I ⇌ FI), so the rate law cannot weight one against the other."
    faces = EnzymeRates.Derivation._faces(am.catalytic)
    built = EnzymeRates.Derivation._built_reading(am)
    eps = EnzymeRates.Derivation._eps_network(
        am, EnzymeRates.Derivation._re_flags(am.catalytic), built)
    w = EnzymeRates.Derivation._admission_face(am, eps,
                                               EnzymeRates.Derivation._law_network(am),
                                               faces, 2)
    @test (w.fault, w.face, w.conformation) == (EnzymeRates.Derivation._PUMPED, 2, 2)
    @test EnzymeRates.Derivation._admission_message(am, w, built) ==
          "Mechanism: at zero P, F + P <--> ES feeds the forms F, FI of " *
          "the inactive conformation, which leave only through the vanishing steps " *
          "F <--> E; FI ⇌ EI while nothing runs back the way it came, so admission " *
          "cannot show that their weight does not depend on how fast those steps vanish."
    # With EI <--> FI steady-state the first built reading stops it, the interior
    # passes, and the face condition rejects the pumped basin at zero P, the ε-path
    # oracle agreeing. The search admits it under a reading in which the :OnlyA steps
    # lose their rates into F's group instead, both oracles agreeing on every face.
    stranded_ss = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + I ⇌ EI        :: NonequalAI
            E + S ⇌ ES        :: OnlyA
            EI <--> FI      :: OnlyA
            F <--> E            :: OnlyA
            F + P <--> ES     :: NonequalAI
            F + I ⇌ FI        :: OnlyA
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
          "Mechanism: at zero P, F + P <--> ES feeds the forms F, FI of " *
          "the inactive conformation, which leave only through the vanishing steps " *
          "F <--> E; FI <--> EI while nothing runs back the way it came, so " *
          "admission cannot show that their weight does not depend on how fast those " *
          "steps vanish."
    @test EnzymeRates.Derivation._guard(stranded_ss) == EnzymeRates.Derivation._WELL_POSED
    @test all(_testhelper_admission_agreement(stranded_ss; rng, corners = true))
    _testhelper_check_law_against_limit(stranded_ss)
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
        _testhelper_check_law_against_limit(am)
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
    _testhelper_check_law_against_limit(bs2)
    # With F <--> E steady-state (R2's graph) the first built reading stops it, the
    # interior passes, and the face condition rejects the pumped basin at zero P; the
    # search admits it as it admits R2, both oracles agreeing on every face.
    bs2_ss = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES        :: OnlyA
            ES <--> EP      :: EqualAI
            F + P ⇌ EP        :: EqualAI
            F <--> E            :: OnlyA
        end
    end
    v = EnzymeRates.Derivation._admission(bs2_ss,
                                          EnzymeRates.Derivation._built_reading(bs2_ss))
    @test (v.fault, v.face, v.conformation) == (PU, 2, 2)
    @test _testhelper_admission_faults(bs2_ss) == [NO, PU, NO]
    @test EnzymeRates.Derivation._guard(bs2_ss) == EnzymeRates.Derivation._WELL_POSED
    @test all(_testhelper_admission_agreement(bs2_ss; rng, corners = true))
    _testhelper_check_law_against_limit(bs2_ss)
end

@testset "a tagged binding that bridges two basins under a reading" begin
    rng = Random.MersenneTwister(59)
    NO, BR = EnzymeRates.Derivation._NO_FAULT, EnzymeRates.Derivation._BRIDGE
    SP, TR = EnzymeRates.Derivation._SPLIT, EnzymeRates.Derivation._TRAP
    # ADbridge: the :OnlyA bindings of I join E and ES of the inactive conformation, which
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
    @test EnzymeRates.Derivation._guard_message(am) ==
          "Mechanism: the :OnlyA step $line bridges forms of the inactive " *
          "conformation that its other steps keep apart, so admission cannot show that " *
          "the limit does not depend on how fast that step vanishes. Make the kinetic " *
          "group of $line steady-state."
    fixed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S <--> ES                         :: EqualAI
            ES <--> E + P                         :: EqualAI
            (E + I <--> EI, ES + I <--> EIS) :: OnlyA
            EI + S ⇌ EIS                      :: EqualAI
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
        row = first(x.catalytic.substrates)
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
    @test EnzymeRates.Derivation._guard_message(drain) ==
          "Mechanism: the :OnlyA step F + I ⇌ FI bridges forms of the " *
          "inactive conformation that its other steps keep apart, so admission cannot " *
          "show that the limit does not depend on how fast that step vanishes. Make the " *
          "kinetic group of F + I ⇌ FI steady-state."
    drain_fixed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES        :: EqualAI
            ES <--> F + P     :: EqualAI
            F <--> E            :: OnlyA
            E + I ⇌ EI        :: OnlyA
            F + I <--> FI     :: OnlyA
            FI ⇌ EI         :: EqualAI
        end
    end
    @test EnzymeRates.Derivation._guard(drain_fixed) == EnzymeRates.Derivation._WELL_POSED
    # R2funnel: under the guard's first built reading, which stops F <--> E, at P zero the
    # trap {F, FS} drains only through a slow step out of a fast transient, and the
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
        _testhelper_check_law_against_limit(x)
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
    row = first(m.substrates)
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
    # R-B1, an equilibrium bridge: the bridge joins E and the dead end EJ, which deletion
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
    @test endswith(EnzymeRates.Derivation._guard_message(am),
                   "Make the kinetic group of " *
                   "$(EnzymeRates.Mechanisms.step_text(am.catalytic, k)) steady-state.")
    fixed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I, J
        steps: begin
            E + S <--> ES         :: EqualAI
            ES <--> E + P         :: EqualAI
            E + J <--> EJ         :: EqualAI
            E + I ⇌ EI            :: OnlyA
            EJ + I <--> EIJ   :: OnlyA
            EI + J ⇌ EIJ      :: EqualAI
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
          "Mechanism: at zero I and P, F + P ⇌ FP runs one way among the " *
          "forms F, FP, FS of the inactive conformation, which leave only through " *
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
        _testhelper_check_law_against_limit(x)
    end
    # R-B3, the set-aside rule's equilibrium condition. In this ping-pong's inactive
    # conformation the cut-off forms EA, FP, F and FB (F is E_res_plusA_minusP)
    # touch free E and EQ, which only the steady-state E + Q <--> EQ joins, so they are
    # not set aside and the interior stays a bridge. Yet the route through them crosses the
    # slow chemistry EA <--> FP, and every corner direction converges there: a
    # conservative rejection.
    D = EnzymeRates.Derivation
    rb3 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                            :: OnlyA
            E + Q <--> EQ                                         :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP      :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP      :: OnlyA
            EA <--> EP_res_plusA_minusP                        :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                           :: OnlyA
        end
    end
    r = fill(D._RUNS, length(rb3.catalytic.steps), 2)
    r[_testhelper_edge(rb3, "E + A ⇌ EA"), 2] = D._LOSES_FORWARD
    r[_testhelper_edge(rb3, "EA <--> EP_res_plusA_minusP"), 2] = D._LOSES_FORWARD
    r[_testhelper_edge(rb3, "E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP"), 2] =
        D._LOSES_REVERSE
    r[_testhelper_edge(rb3, "EB_res_plusA_minusP ⇌ EQ"), 2] = D._LOSES_REVERSE
    # The bridge is led by the last of the two vanishing steps on a cheapest tree, in
    # step order.
    @test D._admission(rb3, r) ==
          D._Verdict(BR, D._LIMIT_UNTESTED, 1, 2,
                     _testhelper_edge(rb3, "E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP"))
    @test _testhelper_epath_verdict(rb3; rng, corners = true,
                                    reading = _testhelper_reading(r)) === :converges
    # R-B4, the set-aside rule's safety clause. In this random-order mechanism's inactive
    # conformation, with P and Q zero, the stopped E + P <--> EP runs one way from EP
    # into free E along its vanishing direction, a pair at a law form that is not safe;
    # the cut-off forms touch free E alone and their fast classes drain, so that pair
    # alone keeps them, and the face stays a bridge. Every corner direction converges
    # there: a conservative rejection.
    rb4 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EP + A ⇌ EAP)        :: OnlyA
            (E + P <--> EP, EA + P <--> EAP)  :: OnlyA
            EA + B <--> EPQ                      :: EqualAI
            EP + Q ⇌ EPQ                         :: OnlyA
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
    @test EnzymeRates.Derivation._admission_message(
        s1, v, EnzymeRates.Derivation._built_reading(s1)) ==
          "Mechanism: at zero A and B, the forms F (inactive) hold the enzyme " *
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
    ad7b_ss = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A
        products: P
        modifiers: I
        steps: begin
            E + A <--> EA                 :: EqualAI
            (E + I ⇌ EI, F + I ⇌ EI)    :: EqualAI
            E + P <--> EP                 :: EqualAI
            E <--> F                        :: OnlyA
            EA <--> EP                  :: EqualAI
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
    # EP leaves only by binding Q, so every other form's weight holds Q, and EP is fed
    # only from EPQ, which EQ + P ⇌ EPQ ties to P: each term of the denominator
    # holds P or Q, so the law is 0/0 with both zero (face 2), where the law's segment
    # chain has two closed classes.
    m = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P, Q
        steps: begin
            E + S ⇌ ES
            ES <--> EPQ
            EQ + P ⇌ EPQ
            E + Q <--> EQ
            EP + Q <--> EPQ
        end
    end
    @test !_testhelper_is_valid(m)
    v = EnzymeRates.Derivation._guard(m)
    @test (v.fault, v.face) == (EnzymeRates.Derivation._LAW_MULTI, 2)
    # The negative control of `_testhelper_finite_on_every_face`, which cannot derive
    # the compiled law of a mechanism the guard rejects: at the concentrations the check
    # reads, the law read without the guard is NaN exactly where P and Q are zero.
    params = Dict{Symbol, Float64}(k => 1.0
                                   for k in EnzymeRates.Derivation._solve(m).fitted)
    params[:Keq] = 2.0
    faces = _testhelper_face_concentrations(m; rng = Random.MersenneTwister(5))
    nan = [isnan(_testhelper_law_value(m, params, Dict(pairs(c)))) for c in faces]
    @test count(nan) == 1
    @test nan == [c.P == 0 && c.Q == 0 for c in faces]
    # The same dead end with EP + Q <--> EPQ in the active conformation only. The
    # active denominator holds P or Q in every term, as above, and the inactive one is
    # weighted by the active free-enzyme weight, which holds Q in every term: the law is
    # 0/0 with P and Q zero, though the inactive denominator alone keeps a constant term.
    # In the inactive conformation EP joins free E only through that :OnlyA step, which
    # the guard's built readings stop, leaving an island; the consistent reading lets it
    # lose its rate into EP, and class 3 rejects zero P and Q, where the law's segment
    # chain has two closed classes, as for the plain mechanism.
    am = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P, Q
        steps: begin
            E + S ⇌ ES          :: EqualAI
            ES <--> EPQ     :: EqualAI
            EQ + P ⇌ EPQ    :: EqualAI
            E + Q <--> EQ       :: EqualAI
            EP + Q <--> EPQ :: OnlyA
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
    # round of expansion (`candidate_mechanisms`) of the bi-bi and ping-pong reactions.
    bibi = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
    end
    pingpong = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
    end
    mechs = EnzymeRates.Mechanism[]
    for rxn in (bibi, pingpong)
        seeds = _testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(rxn))
        append!(mechs, seeds, _testhelper_children(rxn, seeds))
    end
    sample = mechs[randperm(Random.MersenneTwister(6), length(mechs))[1:300]]
    @test isempty(filter(m -> !_testhelper_finite_on_every_face(
                             m; rng = Random.MersenneTwister(5)), sample))
end

@testset "a redundant steady-state binding traps S1-bound enzyme at zero substrates" begin
    rng = Random.MersenneTwister(73)
    NO, TR, SP = EnzymeRates.Derivation._NO_FAULT, EnzymeRates.Derivation._TRAP,
        EnzymeRates.Derivation._SPLIT
    # With S1 and S2 zero, ES1S2 feeds ES1 by the fast release of S2, and ES1
    # leaves only by its slow release of S1: a trap, whose limit depends on the RE speeds.
    # With S2 zero, the lowest forms are joined only through ES1S2: a split.
    m = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S1, S2
        products: P
        steps: begin
            E + S1 <--> ES1
            E + S2 ⇌ ES2
            ES1 + S2 ⇌ ES1S2
            ES2 + S1 ⇌ ES1S2
            ES1S2 <--> EP
            EP ⇌ E + P
        end
    end
    @test _testhelper_face_faults(m) == [NO, NO, TR, NO, SP]
    @test _testhelper_respeed_verdicts(m; rng, law = :limit) ==
          [:agree, :agree, :pattern, :agree, :pattern]
    feed = _testhelper_edge(m, "ES1 + S2 ⇌ ES1S2")
    v = EnzymeRates.Derivation._guard(m)
    @test (v.fault, v.limit, v.face, v.edge) ===
          (TR, EnzymeRates.Derivation._LIMIT_SPEED, 3, feed)
    @test occursin(EnzymeRates.Mechanisms.step_text(m.catalytic, feed),
                   EnzymeRates.Derivation._guard_message(m))
    # The variant moves the steady-state binding of S1 from free E to ES2, where it
    # stays redundant inside the one RE segment: every face passes, and its law is the
    # exact steady state, before any pin.
    spec = only(s for s in MECHANISM_TEST_SPECS
                if s.name == "Numerator: redundant SS-bind")
    x = spec.mechanism
    _testhelper_check_law_against_oracle(spec.mechanism)
    @test _testhelper_face_faults(x) == fill(NO, 5)
    @test _testhelper_respeed_verdicts(spec.mechanism; rng, law = :package) ==
          fill(:agree, 5)
    @test EnzymeRates.Derivation._guard(x) == EnzymeRates.Derivation._WELL_POSED
    @test _testhelper_finite_on_every_face(spec.mechanism; rng = Random.MersenneTwister(5))
    @test length(parameters(spec.mechanism)) == 6
    # The fix the message prints, the feeding step steady-state, passes as well.
    fed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S1, S2
        products: P
        steps: begin
            E + S1 <--> ES1
            E + S2 ⇌ ES2
            ES1 + S2 <--> ES1S2
            ES2 + S1 ⇌ ES1S2
            ES1S2 <--> EP
            EP ⇌ E + P
        end
    end
    @test EnzymeRates.Derivation._guard(fed) === EnzymeRates.Derivation._WELL_POSED
    @test _testhelper_respeed_verdicts(fed; rng) == fill(:agree, 5)
end

@testset "the LDH i-state 6-group mechanism traps active lactate-bound enzyme" begin
    rng = Random.MersenneTwister(79)
    NO, TR, SP = EnzymeRates.Derivation._NO_FAULT, EnzymeRates.Derivation._TRAP,
        EnzymeRates.Derivation._SPLIT
    # With Lactate and NAD zero, ELactateNAD feeds ELactate by the fast release of
    # NAD, and ELactate leaves only slowly: a trap. With NAD zero, a split.
    am = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: NADH, Pyruvate
        products: Lactate, NAD
        oligomeric_state: 4
        steps: begin
            (E + Lactate <--> ELactate,
             ENADH + Lactate <--> ELactateNADH) :: NonequalAI
            (E + NAD ⇌ ENAD,
             ELactate + NAD ⇌ ELactateNAD,
             EPyruvate + NAD ⇌ ENADPyruvate) :: EqualAI
            (E + NADH <--> ENADH,
             ELactate + NADH <--> ELactateNADH,
             EPyruvate + NADH <--> ENADHPyruvate) :: OnlyA
            (E + Pyruvate ⇌ EPyruvate,
             ENAD + Pyruvate ⇌ ENADPyruvate,
             ENADH + Pyruvate ⇌ ENADHPyruvate) :: EqualAI
            ENAD + Lactate ⇌ ELactateNAD :: EqualAI
            ENADHPyruvate <--> ELactateNAD :: OnlyA
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
                   EnzymeRates.Derivation._guard_message(am))
    # The variant makes the :EqualAI NAD group steady-state, keeping its six groups: every
    # face passes the guard and agrees with both oracles, and its law is exact.
    spec = only(s for s in MECHANISM_TEST_SPECS
                if s.name == "LDH i-state NonequalAI 6-group, steady-state NAD group")
    x = spec.mechanism
    _testhelper_check_law_against_oracle(spec.mechanism)
    @test _testhelper_admission_faults(x) == fill(NO, 7)
    @test all(_testhelper_admission_agreement(x; rng))
    @test _testhelper_respeed_verdicts(spec.mechanism; rng, law = :package) ==
          fill(:agree, 7)
    @test EnzymeRates.Derivation._guard(x) == EnzymeRates.Derivation._WELL_POSED
    @test _testhelper_finite_on_every_face(spec.mechanism; rng = Random.MersenneTwister(5))
    @test length(parameters(spec.mechanism)) == 10
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
            E + S ⇌ ES
            ES <--> EP
            EP <--> F + P
        end
    end
    re_cycle = @_testhelper_unguarded @enzyme_mechanism begin
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
    trap = @_testhelper_unguarded @enzyme_mechanism begin
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
    no_limit = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES      :: OnlyA
            ES <--> EP    :: EqualAI
            E + P ⇌ EP      :: EqualAI
        end
    end
    faults = [EnzymeRates.Derivation._guard(m).fault
              for m in (no_cycle, re_cycle, trap, no_limit)]
    @test faults == [EnzymeRates.Derivation._NO_CYCLE, EnzymeRates.Derivation._RE_TURNOVER,
                     EnzymeRates.Derivation._TRAP, EnzymeRates.Derivation._FAST_STEP_PINNED]
    thrown(m) = ErrorException(
        EnzymeRates.Derivation._guard_message(m))
    @test_throws thrown(no_cycle) @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP <--> F + P
        end
    end
    @test_throws thrown(re_cycle) @enzyme_mechanism begin
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
    @test_throws thrown(trap) @enzyme_mechanism begin
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
    @test_throws thrown(no_limit) @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES      :: OnlyA
            ES <--> EP    :: EqualAI
            E + P ⇌ EP      :: EqualAI
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
    both = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                          :: EqualAI
            EA <--> EP_res_plusA_minusP                      :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP    :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP    :: EqualAI
            EB_res_plusA_minusP ⇌ EQ                         :: OnlyA
            E + Q ⇌ EQ                                          :: EqualAI
        end
    end
    @test EnzymeRates.Derivation._guard(both) == EnzymeRates.Derivation._WELL_POSED
    @test EnzymeRates.DSL.parse_mechanism(sprint(show, both)) == both
    no_limit = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EP_res_plusA_minusP                       :: OnlyA
            E + Q <--> EQ                                         :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP      :: EqualAI
            E_res_plusA_minusP + B ⇌ EQ                        :: OnlyA
        end
    end
    v = EnzymeRates.Derivation._guard(no_limit)
    half2 = _testhelper_edge(no_limit, "E_res_plusA_minusP + B ⇌ EQ")
    @test (v.fault, v.edge) == (EnzymeRates.Derivation._FILLED, half2)
    msg = EnzymeRates.Derivation._guard_message(no_limit)
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
                _testhelper_tag_sets[name]; seed = i)
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
    onlya = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES            :: EqualAI
            ES <--> E + P         :: EqualAI
            E + I ⇌ EI            :: OnlyA
        end
    end
    onlyi = fx.competitive_RE
    @test EnzymeRates.Mechanisms.steps_in_active_conformation(onlyi) ==
          EnzymeRates.Mechanisms.steps_in_inactive_conformation(onlya)
    @test EnzymeRates.Mechanisms.steps_in_inactive_conformation(onlyi) ==
          EnzymeRates.Mechanisms.steps_in_active_conformation(onlya)
    @test !all(EnzymeRates.Mechanisms.steps_in_active_conformation(onlyi))
end

@testset "an RE binding whose affinity vanishes still closes a fast cycle" begin
    # The B1 set: an :OnlyI binding at rapid equilibrium loses its affinity in the active
    # conformation, yet on the ε-network it stays fast and closes, with the other
    # rapid-equilibrium steps, a cycle that turns the reaction over: class 2 rejects on
    # the interior. The law networks hold no such cycle, so class 3 and the RE-speed
    # oracle agree on every face. The ε-path finds no path (checked on B1a): the inactive
    # conformation cannot switch off B1a's chemistry of the inhibitor complex, an :OnlyA
    # step at rapid equilibrium.
    b1a = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES            :: EqualAI
            ES <--> EP          :: EqualAI
            E + P ⇌ EP            :: EqualAI
            E + I ⇌ EI            :: OnlyI
            EI + S ⇌ EIS      :: EqualAI
            EIS ⇌ EIP       :: OnlyA
            EP + I ⇌ EIP      :: OnlyI
        end
    end
    b1b = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES            :: OnlyA
            E + P ⇌ EP            :: OnlyA
            ES <--> EP          :: EqualAI
            ES + I ⇌ EIS      :: OnlyI
            EIS ⇌ EIP       :: EqualAI
            EP + I ⇌ EIP      :: OnlyI
        end
    end
    b1d = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES                          :: EqualAI
            ES <--> EP                        :: EqualAI
            E + P ⇌ EP                          :: EqualAI
            (E + I ⇌ EI, EP + I ⇌ EIP)    :: OnlyI
            EI + S ⇌ EIS                    :: EqualAI
            EIS ⇌ EIP                     :: OnlyA
        end
    end
    b1e = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES                    :: EqualAI
            ES <--> EP                  :: EqualAI
            E + P ⇌ EP                    :: EqualAI
            E + P ⇌ EPmod          :: OnlyI
            EPmod + S ⇌ EPmodS    :: EqualAI
            EPmodS ⇌ EPPmod     :: OnlyA
            EP + P ⇌ EPPmod    :: OnlyI
        end
    end
    b1f = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            E + A ⇌ EA                    :: EqualAI
            EA + B ⇌ EAB              :: EqualAI
            EAB <--> EPQ            :: EqualAI
            EQ + P ⇌ EPQ              :: EqualAI
            E + Q ⇌ EQ                    :: EqualAI
            EAB + I ⇌ EABI        :: OnlyI
            EABI ⇌ EIPQ         :: OnlyA
            EPQ + I ⇌ EIPQ        :: OnlyI
        end
    end
    bridge_REcycle = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES                      :: EqualAI
            ES ⇌ F + P                      :: EqualAI
            F <--> E                          :: EqualAI
            (E + I ⇌ EI, F + I ⇌ FI)      :: OnlyI
            FI ⇌ EI                       :: OnlyA
        end
    end
    bridge_REcycle2 = @_testhelper_unguarded @enzyme_mechanism begin
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
    for (i, (label, m)) in enumerate(pairs((; b1a, b1b, b1d, b1e, b1f, bridge_REcycle,
                                               bridge_REcycle2)))
        v = EnzymeRates.Derivation._guard(m)
        @test (label, v.fault, v.face) == (label, EnzymeRates.Derivation._RE_TURNOVER, 0)
        msg = EnzymeRates.Derivation._guard_message(m)
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
    b1c = @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES            :: EqualAI
            ES <--> EP          :: EqualAI
            E + P ⇌ EP            :: EqualAI
            E + I ⇌ EI            :: OnlyI
            EI + S ⇌ EIS      :: EqualAI
            EIS <--> EIP    :: OnlyA
            EP + I ⇌ EIP      :: OnlyI
        end
    end
    @test all(_testhelper_admission_agreement(b1c; rng = Random.MersenneTwister(63)))
    # With bridge_REcycle's I bindings at steady state no fast cycle is left; class 3
    # then rejects zero products, where the steady state depends on how fast the
    # rapid-equilibrium steps run, and passes the other faces. The guard rejects the
    # interior first: the inactive conformation cannot switch off FI ⇌ EI.
    ssbind = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES                          :: EqualAI
            ES ⇌ F + P                          :: EqualAI
            F <--> E                              :: EqualAI
            (E + I <--> EI, F + I <--> FI)    :: OnlyI
            FI ⇌ EI                           :: OnlyA
        end
    end
    v = EnzymeRates.Derivation._guard(ssbind)
    @test (v.fault, v.conformation, v.edge) ==
          (EnzymeRates.Derivation._FAST_STEP_PINNED, 2, _testhelper_edge(ssbind, "FI ⇌ EI"))
    w = _testhelper_face_verdicts(ssbind)[2]
    @test w.face == 2 && w.limit === EnzymeRates.Derivation._LIMIT_SPEED
    @test _testhelper_respeed_verdicts(ssbind; rng = Random.MersenneTwister(64),
                                       law = :limit) ==
          [:agree, :pattern, :agree, :agree, :agree, :agree]
    # With FI <--> EI steady-state as well it is stopped, the interior passes, and the
    # guard itself rejects zero products, where the steady state depends on the speeds;
    # both oracles agree on every face.
    ssbind_ss = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES                          :: EqualAI
            ES ⇌ F + P                          :: EqualAI
            F <--> E                              :: EqualAI
            (E + I <--> EI, F + I <--> FI)    :: OnlyI
            FI <--> EI                        :: OnlyA
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
    # off; but they join the active conformation's group of EIJS and EIJP to
    # free E's, and consistent readings let them lose their rates into it. With the
    # product binding :OnlyA, the inactive conformation keeps the I, J steps, which still
    # join EP to free E: the product binding lies inside one group, the fault no reading
    # changes. No fast cycle is left to class 2.
    m = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I, J
        steps: begin
            E + S ⇌ ES                    :: EqualAI
            ES + I + J ⇌ EIJS       :: OnlyI
            EIJS ⇌ EIJP         :: EqualAI
            EP + I + J ⇌ EIJP       :: OnlyI
            E + P ⇌ EP                    :: OnlyA
        end
    end
    k = _testhelper_edge(m, "EP + I + J ⇌ EIJP")
    @test EnzymeRates.Derivation._admission_interior(
        m, EnzymeRates.Derivation._built_reading(m)) ==
          EnzymeRates.Derivation._Verdict(EnzymeRates.Derivation._FAST_STEP_PINNED,
                                          EnzymeRates.Derivation._LIMIT_UNTESTED, 0, 1, k)
    p = _testhelper_edge(m, "E + P ⇌ EP")
    v = EnzymeRates.Derivation._guard(m)
    @test (v.fault, v.face, v.conformation, v.edge) ==
          (EnzymeRates.Derivation._FAST_STEP_PINNED, 0, 2, p)
    @test EnzymeRates.Derivation._guard_message(m) == _testhelper_pinned_cause(m, p, 2)
    # The same with the chemistry steady-state.
    iso_ss = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I, J
        steps: begin
            E + S ⇌ ES                    :: EqualAI
            ES + I + J ⇌ EIJS       :: OnlyI
            EIJS <--> EIJP      :: EqualAI
            EP + I + J ⇌ EIJP       :: OnlyI
            E + P ⇌ EP                    :: OnlyA
        end
    end
    @test EnzymeRates.Derivation._guard(iso_ss).fault ===
          EnzymeRates.Derivation._FAST_STEP_PINNED
    # With the product binding untagged, the inactive conformation keeps every step: the
    # only cycle of rapid-equilibrium steps, which turns the reaction over, runs there, and
    # class 2, which runs before admission, names that conformation.
    kept = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I, J
        steps: begin
            E + S ⇌ ES                    :: EqualAI
            ES + I + J ⇌ EIJS       :: OnlyI
            EIJS ⇌ EIJP         :: EqualAI
            EP + I + J ⇌ EIJP       :: OnlyI
            E + P ⇌ EP                    :: EqualAI
        end
    end
    v = EnzymeRates.Derivation._guard(kept)
    @test (v.fault, v.face, v.conformation) == (EnzymeRates.Derivation._RE_TURNOVER, 0, 2)
    @test EnzymeRates.Derivation._guard_message(kept) ==
          "Mechanism: the rapid-equilibrium steps EP + I + J ⇌ EIJP, " *
          "E + P ⇌ EP, E + S ⇌ ES, ES + I + J ⇌ EIJS and " *
          "EIJS ⇌ EIJP close a cycle in the inactive conformation that turns " *
          "the reaction over, so the mechanism has no finite rate. Make one step of the " *
          "cycle steady-state."
end

@testset "the ε-network's second closed class rejects the face" begin
    # fz101_1964: at zero products the :OnlyA steps vanish in the inactive conformation
    # and the :OnlyI binding in the active one, and the alive transitions of the
    # ε-network leave two closed classes, so the steady state depends on where the enzyme
    # starts; the same at zero products and modifiers. Class 3 passes every face: the law
    # networks keep one class each.
    fz = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES            :: OnlyA
            ES <--> E + P         :: OnlyA
            ES + I ⇌ EIS      :: OnlyI
            EIS ⇌ EI + P      :: NonequalAI
        end
    end
    v = EnzymeRates.Derivation._guard(fz)
    @test (v.fault, v.face, v.conformation) == (EnzymeRates.Derivation._EPS_NONUNIQUE, 2, 1)
    # The class the message names, EI of the active conformation, which only the
    # vanishing I binding at ES feeds, holds no form of the law networks.
    @test EnzymeRates.Derivation._guard_message(fz) ==
          "Mechanism: at zero P, the forms EI (active) hold the enzyme they " *
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
    # copy_group_bridge and B2b: the :OnlyI bindings of I at E and ES bridge, in the
    # active conformation, forms its other steps keep apart, so the limit depends on how
    # fast the vanishing affinities fall; the ε-path settles on another value. With the
    # bridging group at steady state the guard admits the mechanism, the fix the message
    # prints. At zero modifiers class 3 also rejects copy_group_bridge.
    bridge = _testhelper_onlyi_fixtures().copy_group_bridge
    b2b = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I, J
        steps: begin
            E + S <--> ES                       :: EqualAI
            ES <--> E + P                       :: NonequalAI
            (E + I ⇌ EI, ES + I ⇌ EIS)    :: OnlyI
            EI + S ⇌ EIS                    :: EqualAI
            E + J ⇌ EJ                          :: OnlyA
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
        EnzymeRates.Derivation._guard_message(bridge)
    @test startswith(msg, "Mechanism: the :OnlyI step ") &&
          occursin(" bridges forms of the active conformation that its other steps " *
                   "keep apart, so admission cannot show that the limit does not depend " *
                   "on how fast that step vanishes. Make the kinetic group of ", msg) &&
          endswith(msg, " steady-state.")
    # Along the oracle's five orders the ε-path settles on another value or depends on
    # the order.
    @test _testhelper_epath_verdict(bridge; rng = Random.MersenneTwister(71)) in
          (:wrong_limit, :direction_dependent)
    fixed = @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S <--> ES                           :: EqualAI
            ES <--> E + P                           :: EqualAI
            (E + I <--> EI, ES + I <--> EIS)  :: OnlyI
            EI + S ⇌ EIS                        :: EqualAI
        end
    end
    @test all(_testhelper_admission_agreement(fixed; rng = Random.MersenneTwister(72)))
    # B2a bridges in both conformations, yet admission passes every face; class 3
    # rejects the three faces with zero modifiers, where the limit differs from the law.
    b2a = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I, J
        steps: begin
            E + S <--> ES                       :: EqualAI
            ES <--> E + P                       :: EqualAI
            (E + I ⇌ EI, ES + I ⇌ EIS)    :: OnlyA
            EI + S ⇌ EIS                    :: EqualAI
            (E + J ⇌ EJ, ES + J ⇌ EJS)    :: OnlyI
            EJ + S ⇌ EJS                    :: EqualAI
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
    # X binds only the inactive conformation, and every turnover passes through EX:
    # the active conformation, which must catalyze, turns nothing over. Tagging the X
    # binding :EqualAI or :NonequalAI, the message's fix, passes.
    ea = _testhelper_onlyi_fixtures().essential_activator
    v = EnzymeRates.Derivation._guard(ea)
    @test (v.fault, v.face, v.conformation) ==
          (EnzymeRates.Derivation._NO_ACTIVE_TURNOVER, 0, 1)
    @test EnzymeRates.Derivation._guard_message(ea) ==
          "Mechanism: without its :OnlyI steps E + X ⇌ EX the active " *
          "conformation turns no substrate into product, but it must catalyze the " *
          "reaction. Tag their kinetic groups :EqualAI or :NonequalAI."
    for (i, tag) in enumerate((:EqualAI, :NonequalAI))
        fixed = eval(:(@enzyme_mechanism begin
            substrates: S
            products: P
            modifiers: X
            steps: begin
                E + X ⇌ EX            :: $tag
                EX + S ⇌ ESX      :: EqualAI
                ESX <--> EX + P   :: EqualAI
            end
        end))
        @test all(_testhelper_admission_agreement(fixed;
                                                  rng = Random.MersenneTwister(80 + i)))
    end
end

@testset "an :OnlyI activator whose retags fail is offered no fix" begin
    # X binds E and ES in the inactive conformation only, and P leaves only from
    # ESX, so the active conformation turns nothing over. Tagging the X bindings
    # :EqualAI or :NonequalAI joins E, EX, ESX and ES in one rapid-equilibrium
    # segment that only the X bindings hold together, so each variant splits it at zero
    # X and the message offers no fix.
    am = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: X
        steps: begin
            E + S <--> ES         :: EqualAI
            ES + X ⇌ ESX      :: OnlyI
            ESX ⇌ EX + P      :: EqualAI
            E + X ⇌ EX            :: OnlyI
        end
    end
    v = EnzymeRates.Derivation._guard(am)
    @test (v.fault, v.face, v.conformation) ==
          (EnzymeRates.Derivation._NO_ACTIVE_TURNOVER, 0, 1)
    ks = findall(EnzymeRates.Derivation._excluded_steps(am, 1))
    @test !EnzymeRates.Derivation._retag_passes(am, ks, :EqualAI)
    @test !EnzymeRates.Derivation._retag_passes(am, ks, :NonequalAI)
    @test EnzymeRates.Derivation._guard_message(am) ==
          "Mechanism: without its :OnlyI steps E + X ⇌ EX; ES + X ⇌ " *
          "ESX the active conformation turns no substrate into product, but it " *
          "must catalyze the reaction."
end

@testset "the active conformation's interior faults name the :OnlyI step" begin
    # The mirror of the :OnlyA interior faults, in conformation 1. Each message offers the
    # fixes whose variants pass. The oracles agree with each verdict.
    rng = Random.MersenneTwister(91)
    FP, U = EnzymeRates.Derivation._FAST_STEP_PINNED, EnzymeRates.Derivation._LIMIT_UNTESTED
    # I binds E and ES; without E + I ⇌ EI, the untagged I binding at ES and S
    # binding at EI still join E to EI, so the I binding cannot lose its affinity,
    # and at rapid equilibrium the active conformation cannot switch it off. Neither fix
    # the message offers passes; untagging it would.
    no_limit = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES            :: EqualAI
            ES <--> E + P         :: EqualAI
            E + I ⇌ EI            :: OnlyI
            EI + S ⇌ EIS      :: EqualAI
            ES + I ⇌ EIS      :: EqualAI
        end
    end
    k = _testhelper_edge(no_limit, "E + I ⇌ EI")
    v = EnzymeRates.Derivation._guard(no_limit)
    @test v == EnzymeRates.Derivation._Verdict(FP, U, 0, 1, k)
    @test EnzymeRates.Derivation._guard_message(no_limit) ==
          _testhelper_pinned_cause(no_limit, k, 1)
    @test all(t -> EnzymeRates.Derivation._retag_passes(no_limit, [k], t),
              (:EqualAI, :NonequalAI))
    @test _testhelper_epath_verdict(no_limit; rng) === :no_path
    # I binds E and J binds EI, both :OnlyI; the untagged J binding at E and I binding
    # at EJ join EIJ back to E, so the J binding at EI cannot lose its affinity,
    # and no built reading switches it off. But EI is cut off alone, and the search lets
    # both :OnlyI bindings lose their rates into it, the oracles agreeing.
    square = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I, J
        steps: begin
            E + S ⇌ ES            :: EqualAI
            ES <--> E + P         :: EqualAI
            E + I ⇌ EI            :: OnlyI
            EI + J ⇌ EIJ      :: OnlyI
            E + J ⇌ EJ            :: EqualAI
            EJ + I ⇌ EIJ      :: EqualAI
        end
    end
    k = _testhelper_edge(square, "EI + J ⇌ EIJ")
    @test EnzymeRates.Derivation._admission_interior(
        square, EnzymeRates.Derivation._built_reading(square)) ==
          EnzymeRates.Derivation._Verdict(FP, U, 0, 1, k)
    @test _testhelper_epath_verdict(square; rng) === :no_path
    @test EnzymeRates.Derivation._guard(square) == EnzymeRates.Derivation._WELL_POSED
    @test all(_testhelper_admission_agreement(square; rng, corners = true))
    _testhelper_check_law_against_limit(square)
    # EI is reached only as EIJ gives off J, and the active conformation's other
    # steps make EIJ: no built reading switches that J binding off, though the search
    # lets it lose its rate into EI. But the inactive conformation keeps it, and at zero
    # I and J EI there holds the enzyme it receives: class 3 finds a second closed
    # class, as the RE-speed oracle does.
    j_release = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I, J
        steps: begin
            E + S ⇌ ES            :: EqualAI
            ES <--> E + P         :: EqualAI
            E + J ⇌ EJ            :: EqualAI
            EJ + I ⇌ EIJ      :: EqualAI
            EI + J ⇌ EIJ      :: OnlyI
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
    @test EnzymeRates.Derivation._guard_message(j_release) ==
          "Mechanism: at zero I and J in the inactive conformation, the forms " *
          "{EI} keep enzyme away from free E: every step that leaves them takes up a " *
          "zeroed metabolite (EI + J ⇌ EIJ), so the steady state is not unique."
    # At steady state the built readings stop that J binding, so EI is reached only
    # through a step that vanishes both ways: an island on the :OnlyI binding's unbound
    # side, which collects the active enzyme. The inactive conformation's second closed
    # class at zero I and J rejects it under every reading, as above.
    island = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I, J
        steps: begin
            E + S ⇌ ES            :: EqualAI
            ES <--> E + P         :: EqualAI
            E + J ⇌ EJ            :: EqualAI
            EJ + I ⇌ EIJ      :: EqualAI
            EI + J <--> EIJ   :: OnlyI
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
    @test EnzymeRates.Derivation._guard_message(island) ==
          "Mechanism: at zero I and J in the inactive conformation, the forms " *
          "{EI} keep enzyme away from free E: every step that leaves them takes up a " *
          "zeroed metabolite (EI + J <--> EIJ), so the steady state is not unique."
end

@testset "a pinned step and its mirror get one message, but for conformation and tag" begin
    # Either conformation is offered the same fixes, each with that conformation's tag,
    # when both conformations catalyze in every fix's variant, as here, where the S path
    # stays untagged. I binds E :OnlyI, and the J bindings join EI back to E, so the
    # active conformation cannot switch the I binding off. Tagging :OnlyI the steps that
    # join its sides passes (and so would untagging it, which the message does not
    # offer). The mirror image, the I binding :OnlyA, gets the same message in the
    # inactive conformation, with :OnlyA.
    D = EnzymeRates.Derivation
    swapped(text) = replace(text, "the active conformation" => "the inactive conformation",
                            ":OnlyI" => ":OnlyA")
    square = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I, J
        steps: begin
            E + S ⇌ ES            :: EqualAI
            ES <--> E + P         :: EqualAI
            E + I ⇌ EI            :: OnlyI
            EI + J ⇌ EIJ      :: EqualAI
            E + J ⇌ EJ            :: EqualAI
            EJ + I ⇌ EIJ      :: EqualAI
        end
    end
    k = _testhelper_edge(square, "E + I ⇌ EI")
    joins = [_testhelper_edge(square, t) for t in ("E + J ⇌ EJ", "EJ + I ⇌ EIJ",
                                                   "EI + J ⇌ EIJ")]
    v = D._guard(square)
    @test (v.fault, v.conformation, v.edge) == (D._FAST_STEP_PINNED, 1, k)
    msg = D._guard_message(square)
    @test msg == _testhelper_pinned_cause(square, k, 1) * " Tag :OnlyI the steps that " *
                 "join its sides: E + J ⇌ EJ; EJ + I ⇌ EIJ; EI + J ⇌ EIJ."
    @test all(t -> D._retag_passes(square, [k], t), (:EqualAI, :NonequalAI))
    @test D._retag_passes(square, joins, :OnlyI)
    mirror = _testhelper_mirror(square)
    w = D._guard(mirror)
    @test (w.fault, w.conformation, w.edge) == (D._FAST_STEP_PINNED, 2, k)
    @test D._guard_message(mirror) == swapped(msg)
    # The one rule that tells the states apart: the active conformation must catalyze.
    # The steps that join the sides of this I binding include the uptake of S, which
    # joins E to EP in one step, so tagging them :OnlyI leaves the active conformation no
    # turnover: only the mirror image is offered that fix.
    no_limit = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S <--> EP         :: EqualAI
            EP ⇌ E + P            :: EqualAI
            E + I ⇌ EI            :: OnlyI
            EI + P ⇌ EIP          :: EqualAI
            EP + I ⇌ EIP          :: EqualAI
        end
    end
    k = _testhelper_edge(no_limit, "E + I ⇌ EI")
    e = no_limit.catalytic.steps[k]
    joins = D._kept_path(no_limit, 1, e.from, e.to)
    retagged = copy(no_limit.allosteric.allosteric_states)
    retagged[[no_limit.catalytic.steps[j].kinetic_group for j in joins]] .= :OnlyI
    @test D._guard(EnzymeRates.Mechanism(
              no_limit.catalytic, EnzymeRates.Mechanisms.Allostery(retagged, 1,
              EnzymeRates.Mechanisms.RegulatorySite[]))).fault === D._NO_ACTIVE_TURNOVER
    msg = D._guard_message(no_limit)
    mirror = _testhelper_mirror(no_limit)
    @test D._guard_message(mirror) ==
          swapped(msg) * " Tag :OnlyA the steps that join its sides: " *
          "E + S <--> EP; EP + I ⇌ EIP; EI + P ⇌ EIP."
    # The same graph with the substrate and product roles the other way round. The path
    # the message offers is the first breadth-first in step order (`_kept_path`), here
    # through the S binding E + S ⇌ ES: retagged :OnlyI it is pinned in its turn, and
    # retagged :OnlyA it pins the mirror's, so neither image is offered a fix, though the
    # path through the P release would pass for the mirror (ROADMAP, known issue).
    p_release_path = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES            :: EqualAI
            ES <--> E + P         :: EqualAI
            E + I ⇌ EI            :: OnlyI
            EI + S ⇌ EIS          :: EqualAI
            ES + I ⇌ EIS          :: EqualAI
        end
    end
    k = _testhelper_edge(p_release_path, "E + I ⇌ EI")
    e = p_release_path.catalytic.steps[k]
    @test [EnzymeRates.Mechanisms.step_text(p_release_path.catalytic, j)
           for j in D._kept_path(p_release_path, 1, e.from, e.to)] ==
          ["E + S ⇌ ES", "ES + I ⇌ EIS", "EI + S ⇌ EIS"]
    msg = D._guard_message(p_release_path)
    @test msg == _testhelper_pinned_cause(p_release_path, k, 1)
    @test D._guard_message(_testhelper_mirror(p_release_path)) == swapped(msg)
end

@testset ":OnlyI on a step that does more than bind a modifier" begin
    # Any group may be :OnlyI while the active conformation still catalyzes the reaction.
    # Each mechanism the macro admits has a law equal to the limit on every face, by the
    # ε-path oracle and by the RE-speed oracle against the package's law.
    D = EnzymeRates.Derivation
    rng = Random.MersenneTwister(101)
    # Two isomerization paths from ES to EP; the direct one runs in the inactive
    # conformation alone. Its mirror image, the direct path :OnlyA, is the same model
    # with L inverted, and the guard admits it too.
    two_paths = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES              :: EqualAI
            ES <--> EP            :: OnlyI
            ES ⇌ EstarS           :: EqualAI
            EstarS <--> EstarP    :: EqualAI
            EstarP ⇌ EP           :: EqualAI
            E + P ⇌ EP              :: EqualAI
        end
    end
    @test all(_testhelper_admission_agreement(two_paths; rng))
    @test all(==(:agree), _testhelper_respeed_verdicts(two_paths; rng))
    mirror = _testhelper_mirror(two_paths)
    @test D._guard(mirror) == D._WELL_POSED
    @test all(_testhelper_admission_agreement(mirror; rng))
    @test all(==(:agree), _testhelper_respeed_verdicts(mirror; rng))
    # An exchange step beside a ternary-complex path, the exchange :OnlyI.
    exchange = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA              :: EqualAI
            EA + B <--> EQ + P    :: OnlyI
            EA + B ⇌ EAB        :: EqualAI
            EAB <--> EPQ      :: EqualAI
            EQ + P ⇌ EPQ        :: EqualAI
            E + Q ⇌ EQ              :: EqualAI
        end
    end
    @test all(_testhelper_admission_agreement(exchange; rng))
    @test all(==(:agree), _testhelper_respeed_verdicts(exchange; rng))
    # The only catalytic step :OnlyI: the active conformation cannot catalyze.
    no_turnover = "the active conformation turns no substrate into product"
    @test_throws no_turnover @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES      :: EqualAI
            ES <--> EP    :: OnlyI
            E + P ⇌ EP      :: EqualAI
        end
    end
    # A rapid-equilibrium S binding :OnlyI whose bound form the chemistry and the P
    # binding keep joined to free E: the active conformation cannot switch it off. The
    # message offers no fix: each would leave the active conformation no turnover
    # (untagging the binding passes, but is not one of the fixes it offers).
    pinned_active = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES      :: OnlyI
            ES <--> EP    :: EqualAI
            E + P ⇌ EP      :: EqualAI
        end
    end
    k = _testhelper_edge(pinned_active, "E + S ⇌ ES")
    v = D._guard(pinned_active)
    @test (v.fault, v.conformation, v.edge) == (D._FAST_STEP_PINNED, 1, k)
    msg = D._guard_message(pinned_active)
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
                          for m in _testhelper_mechanisms(
                              EnzymeRates.Enumeration.init_mechanisms(rxn))
                          if !_testhelper_is_allosteric(m)]); by = m -> sprint(show, m))
    rng = Random.MersenneTwister(20261119)
    function _testhelper_random_tagging(holds_only_i::Bool)
        m = rand(rng, plain)
        G = length(_testhelper_group_members(m.catalytic))
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
        EnzymeRates.Mechanism(m.catalytic, EnzymeRates.Mechanisms.Allostery(states, 1,
                                        EnzymeRates.Mechanisms.RegulatorySite[]))
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
                          for m in _testhelper_mechanisms(
                              EnzymeRates.Enumeration.init_mechanisms(rxn))
                          if !_testhelper_is_allosteric(m) &&
                             length(m.catalytic.steps) <= 8]);
                  by = m -> sprint(show, m))
    built(am) = [D._built_reading(am, r) for r in instances(D._SlowUptakeReading)]
    rng = Random.MersenneTwister(20261120)
    sample = EnzymeRates.Mechanism[]
    # The 40 take 1,448 draws; the bound only stops a guard that admits fewer from
    # drawing forever.
    for _ in 1:10_000
        length(sample) == 40 && break
        m = rand(rng, plain)
        G = length(_testhelper_group_members(m.catalytic))
        bits = rand(rng, 1:(2^G - 1))
        am = EnzymeRates.Mechanism(
            m.catalytic,
            EnzymeRates.Mechanisms.Allostery(
                [isodd(bits >> (g - 1)) ? :OnlyA : :EqualAI for g in 1:G], 1,
                EnzymeRates.Mechanisms.RegulatorySite[]))
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
    # conformation reaches its free E without half 2's vanishing steps. M binds EA and
    # the covalent form :OnlyI, and M-bound half 1 runs at rapid equilibrium in the active
    # conformation, where the vanishing M binding at EA lets free E's fast class leave,
    # one way, for the M-bound covalent form. Each conformation's vanishing transitions
    # alone leave the cheapest spanning in-tree at 4 slow transitions; together they cut
    # it to 3. The guard rejects the interior first: the M bindings, which the inactive
    # conformation keeps, join the two sides of M-bound half 1, an :OnlyA step at rapid
    # equilibrium, so that conformation cannot switch it off.
    am = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: M
        steps: begin
            E + A ⇌ EA                                             :: EqualAI
            EA <--> E_res_plusA_minusP + P                      :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP       :: EqualAI
            EB_res_plusA_minusP <--> E + Q                        :: OnlyA
            EA + M ⇌ EAM                                       :: OnlyI
            E_res_plusA_minusP + M <--> EM_res_plusA_minusP    :: OnlyI
            EAM ⇌ EM_res_plusA_minusP + P                     :: OnlyA
        end
    end
    NO, BR, U = EnzymeRates.Derivation._NO_FAULT, EnzymeRates.Derivation._BRIDGE,
        EnzymeRates.Derivation._LIMIT_UNTESTED
    k = _testhelper_edge(am, "EM_res_plusA_minusP + P ⇌ EAM")
    v = EnzymeRates.Derivation._guard(am)
    @test v == EnzymeRates.Derivation._Verdict(EnzymeRates.Derivation._FAST_STEP_PINNED, U,
                                               0, 2, k)
    @test EnzymeRates.Derivation._guard_message(am) ==
          _testhelper_pinned_cause(am, k, 2) * " Make its kinetic group steady-state."
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
            first(am.catalytic.substrates)) ==
              EnzymeRates.Derivation._Verdict(BR, U, 0, 0, 0)
    end
    # The face condition's message names the vanishing steps a cheapest in-tree needs: the
    # M binding at EA of the active conformation and half 2 of the inactive one.
    w = EnzymeRates.Derivation._admission_face(am, eps, law, faces, 2)
    @test w == EnzymeRates.Derivation._Verdict(BR, U, 2, 0, 0)
    @test EnzymeRates.Derivation._admission_message(
        am, w, EnzymeRates.Derivation._built_reading(am)) ==
          "Mechanism: at zero P and Q, the tagged steps EA + M ⇌ EAM; " *
          "E + Q <--> EB_res_plusA_minusP of both conformations together bridge forms " *
          "that each conformation's other steps keep apart, so admission cannot show " *
          "that the limit does not depend on how fast those steps vanish."
end

@testset "_redundant_copy_groups builds each conformation from its own law" begin
    # An aggregate regression pin over the thirteen :OnlyI fixtures: the copy rule reads
    # the active conformation without its :OnlyI groups and the inactive one without its
    # :OnlyA groups, and it flags none of them.
    for (label, am) in pairs(_testhelper_onlyi_fixtures())
        @test (label, EnzymeRates.Enumeration._redundant_copy_groups(
            _testhelper_candidate(am))) == (label, Int[])
    end
end

@testset "check_valid_rate_equation(Bool, m) returns false on a fault" begin
    # A mechanism whose law, on a face (face 2), has several closed classes in its
    # segment chain (`_LAW_MULTI`) fails the guard; the `Bool` form gives the verdict
    # instead of throwing.
    good = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
        end
    end
    bad = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P, Q
        steps: begin
            E + S ⇌ ES
            ES <--> EPQ
            EQ + P ⇌ EPQ
            E + Q <--> EQ
            EP + Q <--> EPQ
        end
    end
    check = EnzymeRates.Derivation.check_valid_rate_equation
    @test check(good) === true
    @test check(Bool, good) === true
    @test check(Bool, bad) === false
    @test @inferred(check(Bool, good)) === true
    @test @inferred(check(Bool, bad)) === false
    v = EnzymeRates.Derivation._guard(bad)
    @test_throws ErrorException(EnzymeRates.Derivation._guard_message(bad)) check(bad)
end
