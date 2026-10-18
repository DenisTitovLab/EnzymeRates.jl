# ABOUTME: Tests of the well-posedness guard on inline fixtures: its networks, faces and
# ABOUTME: classes, their verdicts and messages, each checked against the test oracles.
using Random

@testset "the guard's faces follow the spec's list and the oracle's" begin
    for spec in MECHANISM_TEST_SPECS
        cm = EnzymeRates.catalytic(_testhelper_lift(spec.mechanism))
        faces = EnzymeRates._faces(cm)
        @test faces[1] == 0
        @test [Set(EnzymeRates._face_names(cm, f)) for f in faces] == _testhelper_faces(cm)
        # No face zeroes a substrate with a product: the guard lists no mixed corner.
        @test !any(f -> (n = EnzymeRates._face_names(cm, f);
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
    @test length.(EnzymeRates._faces.((bibi, inhibited, terter))) == (7, 14, 9)
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
    @test EnzymeRates._faces(copy_P) == UInt64[0, 1, 2]
    @test EnzymeRates._face_names(copy_P, UInt64(1)) == [:P]
end

@testset "a law the guard admits can be 0/0 at a mixed corner, which no fit meets" begin
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
    @test EnzymeRates._guard(EnzymeRates.AllostericMechanism(am)) ===
          EnzymeRates._WELL_POSED
    params = _testhelper_random_independent_params(am; rng = Random.MersenneTwister(1))
    @test isnan(rate_equation(am, (A = 1.0, B = 0.0, D = 1.0, P = 0.0, Q = 0.0), params))
    @test rate_equation(am, (A = 0.0, B = 1.0, D = 1.0, P = 0.0, Q = 0.0), params) == 0.0
    # A fit never meets the corner: a zero rate is refused, and the loss scores a NaN
    # prediction as it scores the 0 a well-posed law predicts there.
    _testhelper_table(rate, A, B) =
        (group = [1, 1], Rate = [1.0, rate], A = [1.0, A], B = [1.0, B], D = [1.0, 1.0],
         P = [0.1, 0.0], Q = [0.1, 0.0])
    @test_throws(ErrorException("Zero rate at row 2: log(0) is undefined"),
                 FittingProblem(am, _testhelper_table(0.0, 1.0, 0.0); Keq = 1.0))
    x = [log(params[p]) for p in EnzymeRates.fitted_params(am)]
    at_nan = FittingProblem(am, _testhelper_table(0.5, 1.0, 0.0); Keq = 1.0)
    at_zero = FittingProblem(am, _testhelper_table(0.5, 0.0, 1.0); Keq = 1.0)
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
    @test EnzymeRates._excludes(r2, bind, 2) && !EnzymeRates._excludes(r2, bind, 1)
    @test EnzymeRates._law_mask(r2, 1) == trues(4)
    # The inactive law deletes both tagged steps, which strands every form but E.
    @test EnzymeRates._law_mask(r2, 2) == falses(4)
    law = EnzymeRates._law_network(r2)
    @test (law.nforms, law.ncopies, length(law.from)) == (4, 2, 10)
    @test law.present == [true, true, true, true, true, false, false, false]
    @test law.edge[9:10] == [0, 0] && law.conformation[9:10] == [0, 0]
    @test (law.from[9], law.to[9]) ==
          (EnzymeRates._node(law, 1, 1), EnzymeRates._node(law, 1, 2))
    @test law.fast[9] && law.fast[10] && law.forward[9] && !law.forward[10]
    @test !any(law.vanishing)
    eps = EnzymeRates._eps_network(r2)
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
    eps_re = EnzymeRates._eps_network(r2, trues(4))
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
    law_mm, eps_mm = EnzymeRates._law_network(mm), EnzymeRates._eps_network(mm)
    @test (law_mm.nforms, law_mm.ncopies, length(law_mm.from)) == (2, 1, 4)
    @test all(f -> getfield(law_mm, f) == getfield(eps_mm, f),
              fieldnames(EnzymeRates._Network))
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
    @test !EnzymeRates._has_turnover_cycle(no_turnover, trues(3))
    v = EnzymeRates._guard(no_turnover)
    @test v === EnzymeRates._Verdict(EnzymeRates._NO_CYCLE, EnzymeRates._LIMIT_UNTESTED,
                                     0, 0, 0)
    @test EnzymeRates._guard_message(no_turnover, v) ==
          "Mechanism: no cycle of steps turns substrates into products, so the rate is " *
          "zero everywhere."
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
    v = EnzymeRates._guard(beside)
    @test (v.fault, v.conformation) == (EnzymeRates._RE_TURNOVER, 1)
    @test EnzymeRates._guard_message(beside, v) ==
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
    allo_message = EnzymeRates._guard_message(beside_allo, EnzymeRates._guard(beside_allo))
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
        net = EnzymeRates._eps_network(m)
        cycle = EnzymeRates._re_turnover_transitions(net)
        @test !isempty(cycle) && all(t -> net.fast[t], cycle)
        v = EnzymeRates._guard(m)
        @test v.fault === EnzymeRates._RE_TURNOVER
        msg = EnzymeRates._guard_message(m, v)
        @test all(t -> occursin(EnzymeRates._step_text(m, net.edge[t]), msg), cycle)
        # The fix is printed only where any one step of the cycle made steady-state cures
        # every fast turnover cycle; on sib some of its steps leave a second fast cycle.
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
    @test all(m -> EnzymeRates._guard(m) === EnzymeRates._WELL_POSED, fixed)
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
    @test isempty(EnzymeRates._re_turnover_transitions(EnzymeRates._eps_network(sib_fixed)))
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
        @test EnzymeRates._guard(m) === EnzymeRates._WELL_POSED
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
        cm = EnzymeRates.catalytic(mech)
        names = (substrates = cm.substrates, products = cm.products,
                 modifiers = EnzymeRates.modifiers(cm),
                 regulators = mech isa EnzymeRates.AllostericMechanism ?
                              EnzymeRates.allosteric_regulators(mech) : Symbol[])
        raw = EnzymeRates._raw(cm)
        for k in eachindex(cm.edges)
            op, lhs, rhs = EnzymeRates._step_parts(
                Meta.parse(EnzymeRates._step_text(cm, k)), "@enzyme_mechanism")
            @test EnzymeRates._raw_step(op, lhs, rhs, names, "@enzyme_mechanism") ==
                  raw.steps[k]
        end
    end
    line(m, text) = EnzymeRates._step_text(m, _testhelper_edge(m, text))
    @test line(pingpong, "E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP") ==
          "E(; residual = A - P) + B ⇌ E(B; residual = A - P)"
    @test line(theorell_chance, "EA + B <--> EQ + P") == "E(A) + B <--> E(Q) + P"
    @test line(merged, "EA + B <--> EPQ") == "E(A) + B <--> E(P, Q)"
    @test line(copies, "ES + Pmod ⇌ EPmodS") == "E(S) + P::Mod ⇌ E(P::Mod, S)"
    @test line(copies, "E + I ⇌ EI") == "E + I ⇌ E(I)"
end
