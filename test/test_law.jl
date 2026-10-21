# ABOUTME: Tests of the rate-law engine (src/engine.jl): spanning trees, the tree count,
# ABOUTME: V×τ, and the plain law against the exact steady state, its form and refusals.

@testset "_spanning_trees enumerates each spanning tree once" begin
    trees(n, pairs) = Set(Set.(EnzymeRates._spanning_trees(n, pairs)))
    @test trees(3, [(1, 2), (2, 3), (1, 3)]) == Set([Set([1, 2]), Set([1, 3]), Set([2, 3])])
    @test length(EnzymeRates._spanning_trees(4, [(1, 2), (2, 3), (3, 4), (1, 4)])) == 4
    @test EnzymeRates._spanning_trees(2, [(1, 2), (1, 2)]) == [[1], [2]]
    @test isempty(EnzymeRates._spanning_trees(3, [(1, 2)]))
    @test EnzymeRates._spanning_trees(1, Tuple{Int, Int}[]) == [Int[]]
end

@testset "_tree_count is the BigInt determinant, saturated at typemax(Int)" begin
    complete(n) = [i == j ? n - 1 : -1 for i in 1:n, j in 1:n]
    @test EnzymeRates._tree_count(complete(5)) == 125
    @test EnzymeRates._tree_count(complete(20)) == typemax(Int)
    @test EnzymeRates._tree_count(zeros(Int, 1, 1)) == 1
    @test EnzymeRates._tree_count([1 -1 0; -1 1 0; 0 0 0]) == 0
    @test EnzymeRates._tree_count([1 0 -1; 0 0 0; -1 0 1]) == 0
    rng = Random.MersenneTwister(31)
    for _ in 1:50
        n = rand(rng, 2:14)
        lap = zeros(Int, n, n)
        for a in 1:n, b in (a + 1):n
            rand(rng) < 0.5 || continue
            k = rand(rng, 1:3)
            lap[a, a] += k; lap[b, b] += k; lap[a, b] -= k; lap[b, a] -= k
        end
        exact = det(BigInt.(lap[2:end, 2:end]))
        @test EnzymeRates._tree_count(lap) ==
              (exact > typemax(Int) ? typemax(Int) : Int(exact))
    end
end

@testset "_eq_complexity is segments times spanning trees" begin
    n_plain = 0
    for spec in MECHANISM_TEST_SPECS
        m = _testhelper_lift(spec.mechanism)
        m isa EnzymeRates.Mechanism || continue
        n_plain += 1
        mask = trues(length(m.edges))
        @test EnzymeRates._eq_complexity(m, mask) ===
              Int(_testhelper_laplacian_tree_count(m, mask))
    end
    @test n_plain == 31
    random_bibi = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S1, S2
        products: P1, P2
        steps: begin
            E + S1 <--> E(S1)
            E + S2 <--> E(S2)
            E(S1) + S2 <--> E(S1, S2)
            E(S2) + S1 <--> E(S1, S2)
            E(S1, S2) <--> E(P1, P2)
            E(P1, P2) <--> E(P1) + P2
            E(P1, P2) <--> E(P2) + P1
            E(P1) <--> E + P1
            E(P2) <--> E + P2
        end
    end)
    @test EnzymeRates._eq_complexity(random_bibi) === 336
    one_segment = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E + P
        end
    end)
    @test EnzymeRates._eq_complexity(one_segment) === 1
end

@testset "the plain law equals the exact steady state on every plain spec" begin
    for (i, spec) in enumerate(MECHANISM_TEST_SPECS)
        m = _testhelper_lift(spec.mechanism)
        m isa EnzymeRates.Mechanism || continue
        @testset "$(spec.name)" begin
            L = EnzymeRates._law(m)
            rng = Random.MersenneTwister(500 + i)
            for _ in 1:3
                point, concs, params = _testhelper_law_inputs(m,
                                                              EnzymeRates._solve(m).fitted,
                                                              rng)
                oracle = _testhelper_steady_state_rate(m, point, concs)
                @test abs(_testhelper_evaluate_law(L, merge(params, concs)) - oracle) <=
                      big(10)^-20 * abs(oracle)
            end
        end
    end
end

@testset "a plain law is H times a positive polynomial in lowest terms" begin
    rng = Random.MersenneTwister(7)
    for spec in MECHANISM_TEST_SPECS
        m = _testhelper_lift(spec.mechanism)
        m isa EnzymeRates.Mechanism || continue
        @testset "$(spec.name)" begin
            L = EnzymeRates._law(m)
            @test L.vars.names == [EnzymeRates.metabolites(m);
                                   EnzymeRates._solve(m).fitted; :Keq]
            env = Dict{Symbol, BigFloat}(x => exp(2 * big(rand(rng)) - 1)
                                         for x in L.vars.names)
            H = prod(env[x] for x in EnzymeRates.substrates(m)) -
                prod(env[x] for x in EnzymeRates.products(m)) / env[:Keq]
            @test isapprox(_testhelper_poly_value(L.H, L.vars, env), H;
                           rtol = big(10)^-60)
            @test !iszero(L.Q) && all(>(0), values(L.Q.terms))
            @test all(>(0), values(L.D.terms))
            @test isempty(EnzymeRates._conc_shift((L.Q, L.D), L.vars))
            @test L.D == sum(L.sigma .* L.trees; init = zero(EnzymeRates._Poly)) * L.shift
        end
    end
end

@testset "a conformation's law on its mask is the law of the steps it keeps" begin
    # An inhibitor that binds only the inactive conformation (:OnlyI). The active mask
    # drops its binding and the form E(I), keeping the Michaelis–Menten steps, whose
    # constants both conformations share (:EqualAI); the inactive mask keeps all three
    # steps, the binding under its inactive constant. Both laws number their variables
    # alike: the mechanism's concentrations, the solve's independents, then Keq.
    am = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)            :: EqualAI
            E(S) <--> E + P         :: EqualAI
            E + I ⇌ E(I)            :: OnlyI
        end
    end)
    mm = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E + P
        end
    end)
    competitive = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E + P
            E + I ⇌ E(I)
        end
    end)
    s = EnzymeRates._solve(am)
    cm, mets = EnzymeRates.catalytic(am), EnzymeRates.metabolites(am)
    active = EnzymeRates._law(cm, s, EnzymeRates._active_mask(am); concentrations = mets)
    inactive = EnzymeRates._law(cm, s, EnzymeRates._inactive_mask(am); inactive = true,
                                concentrations = mets)
    @test active.vars == inactive.vars
    @test active.vars.names == [mets; s.names[s.independents]; :Keq]
    rng = Random.MersenneTwister(71)
    for _ in 1:3
        env = Dict{Symbol, BigFloat}(x => exp(2 * big(rand(rng)) - 1)
                                     for x in active.vars.names)
        env[:K_EI_to_E_I] = env[:K_I_EI_to_E_I]
        @test isapprox(_testhelper_evaluate_law(active, env),
                       _testhelper_evaluate_law(EnzymeRates._law(mm), env);
                       rtol = big(10)^-60)
        @test isapprox(_testhelper_evaluate_law(inactive, env),
                       _testhelper_evaluate_law(EnzymeRates._law(competitive), env);
                       rtol = big(10)^-60)
    end
end

@testset "rewriting or regrouping a mechanism changes no plain law" begin
    # Written backwards: every step the other way round, in the reverse order, which
    # canonical form gives back.
    for spec in MECHANISM_TEST_SPECS
        m = _testhelper_lift(spec.mechanism)
        m isa EnzymeRates.Mechanism || continue
        raw = EnzymeRates._raw(m)
        back = EnzymeRates.Mechanism(EnzymeRates.RawMechanism(raw.substrates, raw.products,
            raw.modifiers, reverse([EnzymeRates.RawStep(s.to, s.from, s.released,
                                                        s.consumed, s.rapid_equilibrium)
                                    for s in raw.steps]), reverse(raw.group)))
        @test EnzymeRates._law(back) == EnzymeRates._law(m)
    end
    # A same-model regrouping: the binding square ties B's two dissociation constants,
    # so tying them as one group as well leaves the model and its fitted names alone.
    folded = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))
            E + B ⇌ E(B)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P)
            E(P) ⇌ E + P
        end
    end)
    tied = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))
            (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
            E(A, B) <--> E(P)
            E(P) ⇌ E + P
        end
    end)
    @test folded != tied
    @test EnzymeRates._solve(folded).fitted == EnzymeRates._solve(tied).fitted
    @test EnzymeRates._law(folded) == EnzymeRates._law(tied)
end

@testset "plain laws satisfy detailed balance" begin
    for (i, spec) in enumerate(MECHANISM_TEST_SPECS)
        m = _testhelper_lift(spec.mechanism)
        m isa EnzymeRates.Mechanism || continue
        L = EnzymeRates._law(m)
        rng = Random.MersenneTwister(600 + i)
        env = Dict{Symbol, BigFloat}(x => exp(4 * big(rand(rng)) - 2)
                                     for x in L.vars.names)
        p1 = first(EnzymeRates.products(m))
        env[p1] = env[:Keq] * prod(env[x] for x in EnzymeRates.substrates(m)) /
                  prod((env[x] for x in EnzymeRates.products(m) if x != p1);
                       init = big(1))
        off = copy(env)
        off[p1] /= 10
        @test abs(_testhelper_evaluate_law(L, env)) <=
              big(10)^-30 * abs(_testhelper_evaluate_law(L, off))
    end
end

@testset "the engine refuses a rapid-equilibrium cycle that turns the reaction over" begin
    # E + S ⇌ E(S) ⇌ E(P) ⇌ E + P closes a rapid-equilibrium cycle that turns S into P
    # beside a steady-state route; the macro rejects it (class 2), and the engine,
    # reading it unchecked, refuses to derive a law for it.
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
    msg = "internal error: a rapid-equilibrium cycle of this mechanism turns the " *
          "reaction over (class 2)"
    @test_throws ErrorException(msg) EnzymeRates._law(beside)
    err = @test_throws ErrorException @enzyme_mechanism begin
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
    @test startswith(err.value.msg,
                     "Mechanism: the rapid-equilibrium steps E + P ⇌ E(P), " *
                     "E + S ⇌ E(S) and E(S) ⇌ E(P) close a cycle")
end

@testset "_law refuses a mask whose steady-state steps leave its segments apart" begin
    # The mask keeps two steady-state steps that share no form, so the segment graph has
    # no spanning tree: V×τ is 0, which the cap lets through, and the engine refuses.
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B <--> E(A, B)
            E(A, B) <--> E(P, Q)
            E(P, Q) <--> E(Q) + P
            E(Q) <--> E + Q
        end
    end)
    mask = BitVector([EnzymeRates._step_text(m, k) in ("E + A <--> E(A)",
                                                       "E(Q) + P <--> E(P, Q)")
                      for k in eachindex(m.edges)])
    @test EnzymeRates._eq_complexity(m, mask) == 0
    msg = "internal error: the steady-state steps of this mechanism " *
          "do not join its rapid-equilibrium segments"
    @test_throws ErrorException(msg) EnzymeRates._law(m, EnzymeRates._solve(m), mask)
end

@testset "_law refuses a law over the term cap" begin
    # The mechanism of "Rate equation too large error" (test_rate_eq_derivation.jl):
    # 11 forms, 16 steps, V×τ about 29,000.
    m_manual = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: R1
        steps: begin
            E + A <--> E(A)
            E(A) <--> F(P)
            F(P) <--> F + P
            F + B <--> F(B)
            F(B) <--> E(Q)
            E(Q) <--> E + Q
            E + R1 <--> E(R1)
            E(A) + R1 <--> E(A, R1)
            F(P) + R1 <--> F(P, R1)
            F + R1 <--> F(R1)
            F(B) + R1 <--> F(B, R1)
            E(R1) + A <--> E(A, R1)
            E(A, R1) <--> F(P, R1)
            F(P, R1) <--> F(R1) + P
            F(R1) + B <--> F(B, R1)
            F(B, R1) <--> E(R1) + Q
        end
    end)
    @test EnzymeRates._eq_complexity(m_manual) > EnzymeRates.MAX_RATE_EQUATION_TERMS
    err = @test_throws EnzymeRates.RateEquationTooLargeError EnzymeRates._law(m_manual)
    @test sprint(showerror, err.value) ==
          "Rate equation for this mechanism has more than 5000 polynomial terms (limit: " *
          "5000). Equations this large take a very long time to compile and are " *
          "unlikely to be practically useful for parameter fitting."
end

@testset "an Int64 overflow in a law is an internal error naming the mechanism" begin
    # Under the term cap no coefficient or exponent leaves Int64. A doctored solve whose
    # dependent reads one constant twice, each time to the power typemax(Int), forces an
    # exponent past it.
    mm = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> E(S)
            E(S) <--> E + P
        end
    end)
    s = EnzymeRates._solve(mm)
    d = only(s.dependents)
    f = first(d.factors).first
    big_power = EnzymeRates._Dependent(d.col, d.keq,
                                       [f => typemax(Int) // 1, f => typemax(Int) // 1],
                                       d.label)
    bad = EnzymeRates._Solve(s.columns, s.names, s.column, s.independents, [big_power],
                             s.fitted)
    err = @test_throws ErrorException EnzymeRates._law(mm, bad, trues(length(mm.edges)))
    @test startswith(err.value.msg, "internal error: the rate law of this mechanism " *
                                    "leaves Int64 arithmetic")
    @test endswith(err.value.msg, sprint(show, mm))
end
