# ABOUTME: Tests of src/Derivation/rate_equation_correctness_guard.jl: the guard's entry,
# ABOUTME: verdicts, networks, faces and faults, classes 1 to 3, the limit test, messages.

# ── Faces, networks, classes 1 to 3, the limit test and messages ──────────────────────

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
    runs = EnzymeRates.Mechanisms.steps_running_in_allosteric_conformations(r2)
    @test runs.active == trues(4)
    # The inactive law deletes both tagged steps, which strands every form but E.
    @test runs.inactive == falses(4)
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
    active = EnzymeRates.Mechanisms.steps_running_in_allosteric_conformations(
        consumes_i_inactive).active
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
    # test_rate_equation_correctness_guard.jl "All-RE catalytic cycle raises": its two
    # mechanisms.
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

# ── Bottomless segments, second closed classes and the null-space helper ──────────────

@testset "the macro rejects a rapid-equilibrium segment with no bottom form" begin
    # Random product release from EPQ at rapid equilibrium, but release from EP
    # and EQ at steady state: the RE segment {EP, EPQ, EQ} has weights
    # P·Q : Kq·P : Kp·Q, all zero at P = Q = 0. The RE approximation carries no
    # parameter for how EPQ splits between EP and EQ there, so the rate is
    # 0/0 at zero products: class 3 rejects that face (`_BOTTOMLESS`), and the RE-speed
    # oracle sees the steady state there move with the release speeds.
    random_release = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P, Q
        steps: begin
            E + S ⇌ ES
            E + P <--> EP
            E + Q <--> EQ
            EP + Q ⇌ EPQ
            EQ + P ⇌ EPQ
            ES <--> EPQ
        end
    end
    v = ER.Derivation._guard(random_release)
    @test (v.fault, v.face) == (ER.Derivation._BOTTOMLESS, 2)
    msg = ER.Derivation._guard_message(random_release)
    @test occursin("has no form whose weight keeps the segment's lowest power of " *
                   "each zeroed metabolite", msg)
    @test msg ==
          "Mechanism: at zero P and Q, the rapid-equilibrium segment " *
          "{EP, EPQ, EQ} has no form whose weight keeps the segment's lowest " *
          "power of each zeroed metabolite, so the forms' shares depend on how those " *
          "concentrations vanish and the rate is 0/0 there. Make the kinetic group " *
          "of EP + Q ⇌ EPQ steady-state."
    # The printed fix passes the guard.
    @test ER.Derivation._guard(@_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P, Q
        steps: begin
            E + S ⇌ ES
            E + P <--> EP
            E + Q <--> EQ
            EP + Q <--> EPQ
            EQ + P ⇌ EPQ
            ES <--> EPQ
        end
    end) === ER.Derivation._WELL_POSED
    @test _testhelper_respeed_verdicts(random_release; rng = Random.MersenneTwister(1),
                                       law = :limit)[2] === :pattern
    @test_throws ErrorException(msg) @enzyme_mechanism begin
        substrates: S
        products: P, Q
        steps: begin
            E + S ⇌ ES
            E + P <--> EP
            E + Q <--> EQ
            EP + Q ⇌ EPQ
            EQ + P ⇌ EPQ
            ES <--> EPQ
        end
    end

    # Mirror on the substrate side: the segment {EA, EB, EAB} is 0/0 at
    # A = B = 0 (face 3).
    mirror = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            E + A <--> EA
            E + B <--> EB
            EA + B ⇌ EAB
            EB + A ⇌ EAB
            E + P ⇌ EP
            EAB <--> EP
        end
    end
    v = ER.Derivation._guard(mirror)
    @test (v.fault, v.face) == (ER.Derivation._BOTTOMLESS, 3)
    @test _testhelper_respeed_verdicts(mirror; rng = Random.MersenneTwister(2),
                                       law = :limit)[3] === :pattern
    @test_throws ErrorException(
        ER.Derivation._guard_message(mirror)) @enzyme_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            E + A <--> EA
            E + B <--> EB
            EA + B ⇌ EAB
            EB + A ⇌ EAB
            E + P ⇌ EP
            EAB <--> EP
        end
    end

    # The same segment with A's modifier copy binding EB at rapid equilibrium:
    # EAmodB joins it, but every form still carries A or B (the copy's
    # concentration is A's), so the segment stays bottomless, whichever order the
    # groups are written in.
    copy_last = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            E + A <--> EA
            E + B <--> EB
            EA + B ⇌ EAB
            EB + A ⇌ EAB
            EB + A ⇌ EAmodB
            E + P ⇌ EP
            EAB <--> EP
        end
    end
    copy_first = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            EB + A ⇌ EAmodB
            EAB <--> EP
            E + P ⇌ EP
            EB + A ⇌ EAB
            EA + B ⇌ EAB
            E + B <--> EB
            E + A <--> EA
        end
    end
    v = ER.Derivation._guard(copy_last)
    @test (v.fault, v.face) == (ER.Derivation._BOTTOMLESS, 3)
    @test ER.Derivation._guard(copy_first) == v
    @test copy_first == copy_last
    @test occursin("segment {EA, EAB, EB, EAmodB} has no form",
                   ER.Derivation._guard_message(copy_last))
    @test_throws ErrorException(
        ER.Derivation._guard_message(copy_last)) @enzyme_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            E + A <--> EA
            E + B <--> EB
            EA + B ⇌ EAB
            EB + A ⇌ EAB
            EB + A ⇌ EAmodB
            E + P ⇌ EP
            EAB <--> EP
        end
    end

    # Ordered release: the segment {EP, EPQ} has weights Kq : Q, so EP is its
    # bottom form and the rate is defined at zero products; every face agrees with
    # the RE-speed oracle.
    m_ordered = @enzyme_mechanism begin
        substrates: S
        products: P, Q
        steps: begin
            E + S ⇌ ES
            E + P <--> EP
            EP + Q ⇌ EPQ
            ES <--> EPQ
        end
    end
    @test all(==(:agree),
              _testhelper_respeed_verdicts(m_ordered; rng = Random.MersenneTwister(3)))

    # Mixed abortive complex: the segment {ES, EP, EPS} has weights
    # K·S : K′·P : S·P, zero only at S = P = 0, a corner outside the guard's faces
    # where no turnover is possible anyway. Accepted.
    m_abortive = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> ES
            E + P <--> EP
            ES + P ⇌ EPS
            EP + S ⇌ EPS
            ES <--> EP
        end
    end
    @test all(==(:agree),
              _testhelper_respeed_verdicts(m_abortive; rng = Random.MersenneTwister(4)))
end

@testset "a free form that holds the enzyme it receives is a second closed class" begin
    # Every step that takes up the substrate A is :OnlyA, and the conformations flip at
    # free E, so the inactive conformation keeps E and EP alone: a second free form
    # F, which binds B and joins E only through F + B <--> EB, holds no inactive
    # enzyme. In the active conformation, with every substrate absent, or A, or B, F
    # holds the enzyme it receives, together with EB while B is present, and never
    # returns it to E: a second closed class, which class 3 rejects (`_NONUNIQUE`) and
    # the RE-speed oracle sees.
    s1 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            E + A <--> EA                   :: OnlyA
            EA + B <--> EAB             :: EqualAI
            EAB <--> EP                 :: OnlyA
            EP <--> E + P                   :: EqualAI
            E + A ⇌ EAmod            :: OnlyA
            EAmod + B ⇌ EAmodB      :: EqualAI
            EB + A ⇌ EAmodB      :: EqualAI
            F + B <--> EB                   :: EqualAI
        end
    end
    @test _testhelper_face_faults(s1)[3:5] == fill(ER.Derivation._NONUNIQUE, 3)
    v = ER.Derivation._guard(s1)
    @test (v.fault, v.face) == (ER.Derivation._NONUNIQUE, 3)
    @test _testhelper_respeed_verdicts(s1; rng = Random.MersenneTwister(12),
                                       law = :limit)[3:5] == fill(:nonunique, 3)
    @test_throws ErrorException(
        ER.Derivation._guard_message(s1)) @enzyme_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            E + A <--> EA                   :: OnlyA
            EA + B <--> EAB             :: EqualAI
            EAB <--> EP                 :: OnlyA
            EP <--> E + P                   :: EqualAI
            E + A ⇌ EAmod            :: OnlyA
            EAmod + B ⇌ EAmodB      :: EqualAI
            EB + A ⇌ EAmodB      :: EqualAI
            F + B <--> EB                   :: EqualAI
        end
    end
    # Two variants keep F and pass every face: F <--> E tagged :OnlyA; or F <--> E
    # untagged with the modifier binding of A untagged, so that no :OnlyA binding has
    # its two sides joined by untagged steps.
    f_onlya = @enzyme_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            E + A <--> EA                   :: OnlyA
            EA + B <--> EAB             :: EqualAI
            EAB <--> EP                 :: OnlyA
            EP <--> E + P                   :: EqualAI
            E + A ⇌ EAmod            :: OnlyA
            EAmod + B ⇌ EAmodB      :: EqualAI
            EB + A ⇌ EAmodB      :: EqualAI
            F + B <--> EB                   :: EqualAI
            F <--> E                          :: OnlyA
        end
    end
    copy_untagged = @enzyme_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            E + A <--> EA                   :: OnlyA
            EA + B <--> EAB             :: EqualAI
            EAB <--> EP                 :: OnlyA
            EP <--> E + P                   :: EqualAI
            E + A ⇌ EAmod            :: EqualAI
            EAmod + B ⇌ EAmodB      :: EqualAI
            EB + A ⇌ EAmodB      :: EqualAI
            F + B <--> EB                   :: EqualAI
            F <--> E                          :: EqualAI
        end
    end
    # Both oracles agree on every face: the RE-speed oracle (class 3) and the ε-path
    # (admission of the :OnlyA tags).
    for (i, variant) in enumerate((f_onlya, copy_untagged))
        @test all(_testhelper_admission_agreement(
            variant; rng = Random.MersenneTwister(12 + i)))
    end
end

@testset "the inactive conformation's segments need a bottom form too" begin
    # The :OnlyA P group holds the fused binding EQ + P ⇌ EstarS beside the plain
    # binding E + P ⇌ EP. In the active conformation that binding joins EP,
    # EPQ and EQ to free E at rapid equilibrium, and E is the segment's bottom
    # form. The inactive conformation lacks it and reaches those forms only through
    # the steady-state chemistry and the steady-state Q binding at E, so its segment
    # {EP, EPQ, EQ} holds a product in every form: class 3 finds it bottomless
    # at zero products. The same untagged steps join the two sides of the tagged
    # rapid-equilibrium binding E + P ⇌ EP, so the inactive conformation cannot
    # switch it off: admission rejects the tag set first (`_FAST_STEP_PINNED`), and the
    # ε-path oracle finds no path.
    am = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P, Q
        steps: begin
            E + S ⇌ ES                               :: EqualAI
            ES <--> EPQ                          :: EqualAI
            EPQ ⇌ EP + Q                         :: EqualAI
            EPQ ⇌ EQ + P                         :: EqualAI
            (E + P ⇌ EP, EQ + P ⇌ EstarS)        :: OnlyA
            E + Q <--> EQ                            :: EqualAI
        end
    end
    @test _testhelper_face_faults(am)[2] === ER.Derivation._BOTTOMLESS
    v = ER.Derivation._guard(am)
    @test v.fault === ER.Derivation._FAST_STEP_PINNED
    @test ER.Mechanisms.step_text(am.catalytic, v.edge) == "E + P ⇌ EP"
    @test _testhelper_epath_verdict(am; rng = Random.MersenneTwister(5)) === :no_path
    @test_throws ErrorException(
        ER.Derivation._guard_message(am)) @enzyme_mechanism begin
        substrates: S
        products: P, Q
        steps: begin
            E + S ⇌ ES                               :: EqualAI
            ES <--> EPQ                          :: EqualAI
            EPQ ⇌ EP + Q                         :: EqualAI
            EPQ ⇌ EQ + P                         :: EqualAI
            (E + P ⇌ EP, EQ + P ⇌ EstarS)        :: OnlyA
            E + Q <--> EQ                            :: EqualAI
        end
    end
end

@testset "the rational nullspace helper" begin
    R = Rational{BigInt}

    @testset "_rational_nullspace" begin
        # Zero map: every x solves M * x = 0, so the nullspace is the whole
        # 3-dimensional space.
        M = zeros(R, 2, 3)
        N = ER.Derivation._rational_nullspace(M)
        @test size(N, 2) == 3

        # Full-rank square matrix: only x = 0 solves M * x = 0.
        M = R[1 2; 3 4]
        N = ER.Derivation._rational_nullspace(M)
        @test size(N, 2) == 0

        # Rank-1 2x2 (row 2 is twice row 1): both rows reduce to x + y = 0,
        # spanned by (1, -1). n = 2, rank = 1.
        M = R[1 1; 2 2]
        N = ER.Derivation._rational_nullspace(M)
        @test M * N == zeros(R, 2, size(N, 2))
        @test size(N, 2) == 2 - 1

        # More rows than columns, rank-deficient: every row is a multiple of
        # (1, 2), i.e. x + 2y = 0, spanned by (2, -1).
        M = R[1 2; 2 4; 3 6]
        N = ER.Derivation._rational_nullspace(M)
        @test M * N == zeros(R, 3, size(N, 2))
        @test size(N, 2) == 2 - 1

        # More columns than rows: x + z = 0 and y + z = 0, spanned by
        # (-1, -1, 1).
        M = R[1 0 1; 0 1 1]
        N = ER.Derivation._rational_nullspace(M)
        @test M * N == zeros(R, 2, size(N, 2))
        @test size(N, 2) == 3 - 2
    end
end

@testset "every derivation entry rejects a bottomless rapid-equilibrium segment" begin
    # Random product release at rapid equilibrium between steady-state releases leaves
    # the segment {EP, EPQ, EQ} without a bottom form, so the rate is 0/0 at zero
    # products. The macro's boundary check (`check_valid_rate_equation`) rejects the
    # mechanism.
    # The constructor checks the representation only, so it builds the mechanism from
    # the walker's raw value, and the mechanism compiles; every derivation entry then
    # rejects it with the same message.
    raw = EnzymeRates.DSL.raw_mechanism(:(begin
        substrates: S
        products: P, Q
        steps: begin
            E + S ⇌ ES
            ES <--> EPQ
            EPQ ⇌ EQ + P
            EPQ ⇌ EP + Q
            EP <--> E + P
            EQ <--> E + Q
        end
    end))
    m = EnzymeRates.DSL._construct(raw)
    v = EnzymeRates.Derivation._guard(m)
    @test (v.fault, v.face) == (EnzymeRates.Derivation._BOTTOMLESS, 2)
    msg = ErrorException(EnzymeRates.Derivation._guard_message(m))
    @test_throws msg EnzymeRates.DSL._checked_mechanism(raw)
    rate = compile_rate_equation(m)
    @test_throws msg parameters(rate)
    @test_throws msg parameters(m)
    @test_throws msg rate((S = 1.0, P = 1.0, Q = 1.0), (;), 1.0)
    @test_throws msg rate_equation_string(m)
    @test_throws msg dependent_parameter_values(m, (;), 1.0)
    @test_throws msg EnzymeRates.Derivation.parameter_units(m)
end

# ── A fast cycle, a trap and shared names ─────────────────────────────────────────────

@testset "the third reproducer with a rapid-equilibrium Lactate release traps ENAD" begin
    # Written with ENAD + Lactate ⇌ ELactateNAD at rapid equilibrium beside the
    # steady-state NAD binding, the mechanism fills ENAD by the fast Lactate release at
    # zero products, and only steady-state steps empty it: a trap (`_TRAP`), which the
    # RE-speed oracle sees as a pattern, as it does the split at zero Lactate. The
    # reproducer takes the fix the message prints, that step steady-state, which passes
    # every face and keeps the steady-state speed k_I_E_NAD_to_ENAD dependent.
    NO = EnzymeRates.Derivation._NO_FAULT
    trap = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: NADH, Pyruvate
        products: Lactate, NAD
        oligomeric_state: 4
        steps: begin
            E + Lactate ⇌ ELactate :: EqualAI
            E + NAD <--> ENAD :: NonequalAI
            E + NADH ⇌ ENADH :: OnlyA
            (E + Pyruvate <--> EPyruvate,
             ENADH + Pyruvate <--> ENADHPyruvate) :: EqualAI
            ELactate + NAD ⇌ ELactateNAD :: EqualAI
            (ELactate + NADH ⇌ ELactateNADH,
             EPyruvate + NADH ⇌ ENADHPyruvate) :: OnlyA
            ENAD + Lactate ⇌ ELactateNAD :: EqualAI
            ENAD + Pyruvate <--> ENADPyruvate :: EqualAI
            ENADH + Lactate ⇌ ELactateNADH :: EqualAI
            ENADHPyruvate <--> ELactateNAD :: OnlyA
            EPyruvate + NAD ⇌ ENADPyruvate :: EqualAI
        end
    end
    v = EnzymeRates.Derivation._guard(trap)
    @test (v.fault, v.limit, v.face) ==
          (EnzymeRates.Derivation._TRAP, EnzymeRates.Derivation._LIMIT_SPEED, 2)
    @test EnzymeRates.Mechanisms.step_text(trap.catalytic, v.edge) ==
          "ENAD + Lactate ⇌ ELactateNAD"
    @test _testhelper_admission_faults(trap) ==
          [NO, EnzymeRates.Derivation._TRAP, NO, NO, NO, EnzymeRates.Derivation._SPLIT, NO]
    @test _testhelper_respeed_verdicts(trap; rng = Random.MersenneTwister(1),
                                       law = :limit) ==
          [:agree, :pattern, :agree, :agree, :agree, :pattern, :agree]
    msg = EnzymeRates.Derivation._guard_message(trap)
    @test_throws ErrorException(msg) @enzyme_mechanism begin
        substrates: NADH, Pyruvate
        products: Lactate, NAD
        oligomeric_state: 4
        steps: begin
            E + Lactate ⇌ ELactate :: EqualAI
            E + NAD <--> ENAD :: NonequalAI
            E + NADH ⇌ ENADH :: OnlyA
            (E + Pyruvate <--> EPyruvate,
             ENADH + Pyruvate <--> ENADHPyruvate) :: EqualAI
            ELactate + NAD ⇌ ELactateNAD :: EqualAI
            (ELactate + NADH ⇌ ELactateNADH,
             EPyruvate + NADH ⇌ ENADHPyruvate) :: OnlyA
            ENAD + Lactate ⇌ ELactateNAD :: EqualAI
            ENAD + Pyruvate <--> ENADPyruvate :: EqualAI
            ENADH + Lactate ⇌ ELactateNADH :: EqualAI
            ENADHPyruvate <--> ELactateNAD :: OnlyA
            EPyruvate + NAD ⇌ ENADPyruvate :: EqualAI
        end
    end
    fixed = ALLOSTERIC_UNDEFVAR_REPRODUCERS[3]
    @test _testhelper_admission_faults(fixed) == fill(NO, 7)
    @test _testhelper_respeed_verdicts(fixed; rng = Random.MersenneTwister(2)) ==
          fill(:agree, 7)
    @test :k_I_E_NAD_to_ENAD in keys(_testhelper_dependent_param_exprs(fixed)[1])
    _testhelper_check_law_against_oracle(ALLOSTERIC_UNDEFVAR_REPRODUCERS[3])
end

@testset "All-RE catalytic cycle raises" begin
    # A cycle of rapid-equilibrium steps that turns the reaction over gives the mechanism
    # no finite rate. The guard's class 2 finds it, and the macro throws when its expanded
    # code runs, naming the cycle's steps.
    msg = "close a cycle that turns the reaction over, so the mechanism has no finite rate"
    # Binding stage mixed (S2 SS, S1 RE), release stage mixed (P1 SS, P2 RE), chemistry
    # RE: a complete all-RE catalytic cycle.
    @test_throws msg @enzyme_mechanism begin
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
    # The activator route ER + S ⇌ ERS ⇌ EPR ⇌ ER + P is an all-RE catalytic
    # cycle.
    @test_throws msg @enzyme_mechanism begin
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
            (E + R ⇌ ER,
             ES + R ⇌ ERS,
             EP + R ⇌ EPR)
        end
    end
end

# Every constant of a rate equation needs a name of its own. An :EqualAI step's constants
# are bare, while a :NonequalAI step's name their conformation, k_A_… for the active one
# and k_I_… for the inactive one. Form names are opaque, so a form named A_ES or I_ES
# renders the same name as ES with a conformation prefix, and an :EqualAI step from it can
# take a :NonequalAI step's name. A metabolite named minusQ imitates a residual. The
# constructors build such mechanisms; the guard rejects them, so the macros do too.

"""The block of an allosteric uni-uni mechanism whose chemistry ES <--> EP is :NonequalAI,
with a second S-bound form named `<prefix>_ES` whose own chemistry to EP is :EqualAI."""
_testhelper_second_bound_block(prefix) = :(begin
    substrates: S
    products:   P
    steps: begin
        E + S ⇌ ES                          :: EqualAI
        ES <--> EP                          :: NonequalAI
        E + S ⇌ $(Symbol(prefix, "_ES"))     :: EqualAI
        $(Symbol(prefix, "_ES")) <--> EP     :: EqualAI
        EP ⇌ E + P                          :: EqualAI
    end
end)

@testset "constants that would share a name are rejected by the guard" begin
    _testhelper_built(prefix) = _testhelper_build(
        EnzymeRates.DSL.raw_mechanism(_testhelper_second_bound_block(prefix)))
    _testhelper_clash(name) = ErrorException(
        "Mechanism: constants of this mechanism's rate equation share a name ($name). " *
        "Rename a metabolite or an enzyme form so that each constant has a name of its " *
        "own.")
    # The :EqualAI step A_ES <--> EP has the constant k_A_ES_to_EP, the name of the
    # active constant of the :NonequalAI step ES <--> EP.
    conformation_a = _testhelper_built(:A)
    clash_a = _testhelper_clash("k_A_ES_to_EP")
    @test !EnzymeRates.Derivation.check_valid_rate_equation(Bool, conformation_a)
    @test_throws clash_a EnzymeRates.Derivation.check_valid_rate_equation(conformation_a)
    @test_throws clash_a eval(:(@enzyme_mechanism $(_testhelper_second_bound_block(:A))))
    @test_throws clash_a parameters(conformation_a)
    @test_throws clash_a rate_equation_string(conformation_a)
    params = (K_ES_to_E_S = 1.0, K_A_ES_to_E_S = 2.0, K_EP_to_E_P = 3.0,
              k_A_ES_to_EP = 4.0, k_I_ES_to_EP = 0.5, L = 10.0)
    @test_throws clash_a dependent_parameter_values(conformation_a, params, 4.0)
    # The compiled law derives on its first call, which raises the guard's error.
    @test_throws clash_a compile_rate_equation(conformation_a)((S = 1.0, P = 0.1),
                                                               params, 4.0)

    # Named I_ES, the :EqualAI constant k_I_ES_to_EP is the name of the :NonequalAI step's
    # inactive constant.
    conformation_i = _testhelper_built(:I)
    clash_i = _testhelper_clash("k_I_ES_to_EP")
    @test !EnzymeRates.Derivation.check_valid_rate_equation(Bool, conformation_i)
    @test_throws clash_i parameters(conformation_i)
    @test_throws clash_i rate_equation_string(conformation_i)

    # Named F_ES, each constant has a name of its own.
    renamed = eval(:(@enzyme_mechanism $(_testhelper_second_bound_block(:F))))
    @test parameters(renamed) == (:K_EP_to_E_P, :K_ES_to_E_S, :K_F_ES_to_E_S, :k_A_ES_to_EP,
                                  :k_F_ES_to_EP, :k_I_ES_to_EP, :L)
    @test renamed == _testhelper_built(:F)

    # A modifier named minusQ imitates a residual: the side E_res_plusA_minusP + minusQ
    # and the form E_res_plusA_minusP_minusQ both render E_res_plusA_minusP_minusQ, so each
    # constant of the minusQ binding takes the name of a constant of the step between
    # EminusQ_res_plusA_minusP and that form; on no cycle, that step runs away from free
    # E, from EminusQ_res_plusA_minusP. The macro runs the guard, so it rejects the
    # mechanism; the constructors build it.
    residual_block = :(begin
        substrates: A, B
        products: P, Q
        modifiers: minusQ
        steps: begin
            E + A ⇌ EA
            EA <--> EP_res_plusA_minusP
            EP_res_plusA_minusP ⇌ E_res_plusA_minusP + P
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
            EB_res_plusA_minusP <--> EQ
            EQ ⇌ E + Q
            E_res_plusA_minusP + minusQ <--> EminusQ_res_plusA_minusP
            EminusQ_res_plusA_minusP <--> E_res_plusA_minusP_minusQ
        end
    end)
    residual_clash = _testhelper_build(EnzymeRates.DSL.raw_mechanism(residual_block))
    clash_residual = ErrorException(
        "Mechanism: constants of this mechanism's rate equation share a name " *
        "(k_E_res_plusA_minusP_minusQ_to_EminusQ_res_plusA_minusP, " *
        "k_EminusQ_res_plusA_minusP_to_E_res_plusA_minusP_minusQ). Rename a metabolite " *
        "or an enzyme form so that each constant has a name of its own.")
    @test !EnzymeRates.Derivation.check_valid_rate_equation(Bool, residual_clash)
    @test_throws clash_residual eval(:(@enzyme_mechanism $residual_block))
    @test_throws clash_residual parameters(residual_clash)
    # The compiled law derives on its first call, which raises the guard's error.
    names = Tuple(metabolites(residual_clash))
    concs = NamedTuple{names}(ntuple(_ -> 1.0, length(names)))
    @test_throws clash_residual compile_rate_equation(residual_clash)(concs, (;), 2.0)
end
