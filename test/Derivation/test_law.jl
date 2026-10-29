# ABOUTME: Tests of the rate-law engine (src/Derivation/engine.jl): spanning trees, V×τ,
# ABOUTME: plain and MWC laws against exact steady states, form and refusals, summed cap.

@testset "_spanning_trees enumerates each spanning tree once" begin
    trees(n, pairs) = Set(Set.(EnzymeRates.Derivation._spanning_trees(n, pairs)))
    @test trees(3, [(1, 2), (2, 3), (1, 3)]) == Set([Set([1, 2]), Set([1, 3]), Set([2, 3])])
    @test length(
        EnzymeRates.Derivation._spanning_trees(4, [(1, 2), (2, 3), (3, 4), (1, 4)])) == 4
    @test EnzymeRates.Derivation._spanning_trees(2, [(1, 2), (1, 2)]) == [[1], [2]]
    @test isempty(EnzymeRates.Derivation._spanning_trees(3, [(1, 2)]))
    @test EnzymeRates.Derivation._spanning_trees(1, Tuple{Int, Int}[]) == [Int[]]
end

@testset "_tree_count is the BigInt determinant, saturated at typemax(Int)" begin
    complete(n) = [i == j ? n - 1 : -1 for i in 1:n, j in 1:n]
    @test EnzymeRates.Derivation._tree_count(complete(5)) == 125
    @test EnzymeRates.Derivation._tree_count(complete(20)) == typemax(Int)
    @test EnzymeRates.Derivation._tree_count(zeros(Int, 1, 1)) == 1
    @test EnzymeRates.Derivation._tree_count([1 -1 0; -1 1 0; 0 0 0]) == 0
    @test EnzymeRates.Derivation._tree_count([1 0 -1; 0 0 0; -1 0 1]) == 0
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
        @test EnzymeRates.Derivation._tree_count(lap) ==
              (exact > typemax(Int) ? typemax(Int) : Int(exact))
    end
end

@testset "rate_equation_complexity is segments times spanning trees" begin
    n_plain = 0
    for spec in MECHANISM_TEST_SPECS
        m = _testhelper_lift(spec.mechanism)
        m isa EnzymeRates.Mechanism || continue
        n_plain += 1
        mask = trues(length(m.steps))
        @test EnzymeRates.Derivation.rate_equation_complexity(m, mask) ===
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
    @test EnzymeRates.Derivation.rate_equation_complexity(random_bibi) === 336
    one_segment = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E + P
        end
    end)
    @test EnzymeRates.Derivation.rate_equation_complexity(one_segment) === 1
end

@testset "the plain law equals the exact steady state on every plain spec" begin
    for (i, spec) in enumerate(MECHANISM_TEST_SPECS)
        m = _testhelper_lift(spec.mechanism)
        m isa EnzymeRates.Mechanism || continue
        @testset "$(spec.name)" begin
            L = EnzymeRates.Derivation._law(m)
            rng = Random.MersenneTwister(500 + i)
            for _ in 1:3
                point, concs, params =
                    _testhelper_law_inputs(m, EnzymeRates.Derivation._solve(m).fitted, rng)
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
            L = EnzymeRates.Derivation._law(m)
            @test L.vars.names == [EnzymeRates.metabolites(m);
                                   EnzymeRates.Derivation._solve(m).fitted; :Keq]
            env = Dict{Symbol, BigFloat}(x => exp(2 * big(rand(rng)) - 1)
                                         for x in L.vars.names)
            H = prod(env[x] for x in m.substrates) -
                prod(env[x] for x in m.products) / env[:Keq]
            @test isapprox(_testhelper_poly_value(L.H, L.vars, env), H;
                           rtol = big(10)^-60)
            @test !iszero(L.Q) && all(>(0), values(L.Q.terms))
            @test all(>(0), values(L.D.terms))
            @test isempty(EnzymeRates.Derivation._conc_shift((L.Q, L.D), L.vars))
            @test L.D == sum(L.sigma .* L.trees;
                             init = zero(EnzymeRates.Derivation._Poly)) * L.shift
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
    s = EnzymeRates.Derivation._solve(am)
    cm, mets = EnzymeRates.Mechanisms.catalytic(am), EnzymeRates.metabolites(am)
    active = EnzymeRates.Derivation._law(
        cm, s, EnzymeRates.Derivation.conformation_steps(am, :active);
        concentrations = mets)
    inactive = EnzymeRates.Derivation._law(
        cm, s, EnzymeRates.Derivation.conformation_steps(am, :inactive);
        inactive = true, concentrations = mets)
    @test active.vars == inactive.vars
    @test active.vars.names == [mets; s.names[s.independents]; :Keq]
    rng = Random.MersenneTwister(71)
    for _ in 1:3
        env = Dict{Symbol, BigFloat}(x => exp(2 * big(rand(rng)) - 1)
                                     for x in active.vars.names)
        env[:K_EI_to_E_I] = env[:K_I_EI_to_E_I]
        @test isapprox(_testhelper_evaluate_law(active, env),
                       _testhelper_evaluate_law(EnzymeRates.Derivation._law(mm), env);
                       rtol = big(10)^-60)
        @test isapprox(
            _testhelper_evaluate_law(inactive, env),
            _testhelper_evaluate_law(EnzymeRates.Derivation._law(competitive), env);
            rtol = big(10)^-60)
    end
end

@testset "rewriting or regrouping a mechanism changes no plain law" begin
    # Written backwards: every step the other way round, in the reverse order, which
    # canonical form gives back.
    for spec in MECHANISM_TEST_SPECS
        m = _testhelper_lift(spec.mechanism)
        m isa EnzymeRates.Mechanism || continue
        raw = EnzymeRates.Mechanisms.RawMechanism(m)
        back = EnzymeRates.Mechanism(EnzymeRates.Mechanisms.RawMechanism(
                                         raw.substrates, raw.products, raw.modifiers,
                                         reverse([EnzymeRates.Mechanisms.RawStep(
                                             s.to, s.from, s.released, s.consumed,
                                             s.rapid_equilibrium) for s in raw.steps]),
                                         reverse(raw.group)))
        @test EnzymeRates.Derivation._law(back) == EnzymeRates.Derivation._law(m)
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
    @test EnzymeRates.Derivation._solve(folded).fitted ==
          EnzymeRates.Derivation._solve(tied).fitted
    @test EnzymeRates.Derivation._law(folded) == EnzymeRates.Derivation._law(tied)
end

@testset "plain laws satisfy detailed balance" begin
    for (i, spec) in enumerate(MECHANISM_TEST_SPECS)
        m = _testhelper_lift(spec.mechanism)
        m isa EnzymeRates.Mechanism || continue
        L = EnzymeRates.Derivation._law(m)
        rng = Random.MersenneTwister(600 + i)
        env = Dict{Symbol, BigFloat}(x => exp(4 * big(rand(rng)) - 2)
                                     for x in L.vars.names)
        p1 = first(m.products)
        env[p1] = env[:Keq] * prod(env[x] for x in m.substrates) /
                  prod((env[x] for x in m.products if x != p1);
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
    @test_throws ErrorException(msg) EnzymeRates.Derivation._law(beside)
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
    mask = BitVector([EnzymeRates.Mechanisms.step_text(m, k) in ("E + A <--> E(A)",
                                                                  "E(Q) + P <--> E(P, Q)")
                      for k in eachindex(m.steps)])
    @test EnzymeRates.Derivation.rate_equation_complexity(m, mask) == 0
    msg = "internal error: the steady-state steps of this mechanism " *
          "do not join its rapid-equilibrium segments"
    @test_throws ErrorException(msg) EnzymeRates.Derivation._law(
        m, EnzymeRates.Derivation._solve(m), mask)
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
    @test EnzymeRates.Derivation.rate_equation_complexity(m_manual) >
          EnzymeRates.Derivation.MAX_RATE_EQUATION_TERMS
    err = @test_throws(EnzymeRates.Derivation.RateEquationTooLargeError,
                       EnzymeRates.Derivation._law(m_manual))
    @test sprint(showerror, err.value) ==
          "Rate equation for this mechanism has more than 5000 polynomial terms (limit: " *
          "5000). Equations this large take a very long time to compile and are " *
          "unlikely to be practically useful for parameter fitting."
end

@testset "an Int64 overflow in a law is an internal error naming the mechanism" begin
    # Under the term cap no coefficient or exponent leaves Int64. A doctored solve whose
    # first dependent reads one constant twice, each time to the power typemax(Int),
    # forces an exponent past it.
    function _testhelper_overflowing(s)
        d = first(s.dependents)
        f = first(d.factors).first
        big_power = EnzymeRates.Derivation._Dependent(
            d.col, d.keq, [f => typemax(Int) // 1, f => typemax(Int) // 1], d.label)
        EnzymeRates.Derivation._Solve(s.columns, s.names, s.column, s.independents,
                                      [big_power; s.dependents[2:end]], s.fitted)
    end
    mm = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> E(S)
            E(S) <--> E + P
        end
    end)
    bad = _testhelper_overflowing(EnzymeRates.Derivation._solve(mm))
    err = @test_throws ErrorException EnzymeRates.Derivation._law(mm, bad,
                                                                  trues(length(mm.steps)))
    @test startswith(err.value.msg, "internal error: the rate law of this mechanism " *
                                    "leaves Int64 arithmetic")
    @test endswith(err.value.msg, sprint(show, mm))
    # An MWC law's message names the allosteric mechanism, not its catalytic graph.
    am = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)       :: EqualAI
            E(S) <--> E(P)     :: NonequalAI
            E + P ⇌ E(P)       :: EqualAI
        end
    end)
    bad = _testhelper_overflowing(EnzymeRates.Derivation._solve(am))
    err = @test_throws ErrorException EnzymeRates.Derivation._law(am, bad)
    @test startswith(err.value.msg, "internal error: the rate law of this mechanism " *
                                    "leaves Int64 arithmetic")
    @test endswith(err.value.msg, sprint(show, am))
end

# ── The MWC law ─────────────────────────────────────────────────────────────

@testset "the MWC law equals the one-subunit oracle on every allosteric spec" begin
    # At one catalytic subunit the MWC law is the exact steady state of one protomer that
    # flips at free E: 3 points per spec, 1e-20 relative in 256-bit BigFloat. Each
    # conformation's Q has positive coefficients, the law vanishes at equilibrium in each
    # conformation and combined, and the fitted names count the declared dimension. A
    # second derivation gives an equal law with an equal hash.
    specs = [s for s in MECHANISM_TEST_SPECS
             if s.mechanism isa EnzymeRates.AllostericEnzymeMechanism]
    @test length(specs) == 11
    for (i, spec) in enumerate(specs)
        @testset "$(spec.name)" begin
            am = _testhelper_lift(_testhelper_single_subunit(spec.mechanism))
            s = EnzymeRates.Derivation._solve(am)
            A = EnzymeRates.Derivation._law(am, s)
            rng = Random.MersenneTwister(i)
            for _ in 1:3
                point, concs, params = _testhelper_law_inputs(am, s.fitted, rng)
                oracle = _testhelper_steady_state_rate(am, point, concs)
                law = _testhelper_evaluate_law(A, merge(params, concs))
                @test abs(law - oracle) <= big(10)^-20 * abs(oracle)
            end
            @test all(>(0), values(A.active.Q.terms))
            @test all(>(0), values(A.inactive.Q.terms))
            @test all(<(big(10)^-30),
                      _testhelper_detailed_balance(am, Random.MersenneTwister(100 + i)))
            @test length(s.fitted) == _testhelper_declared_dimension(am)
            again = EnzymeRates.Derivation._law(am, s)
            @test again == A && hash(again) == hash(A)
        end
    end
end

@testset "the MWC law equals the n-protomer oracle at 1–3 protomers" begin
    # The bi-uni of the :NonequalAI ground-truth gate: A binds at steady state, B and P at
    # rapid equilibrium, and the chemistry differs between the conformations, so the
    # inactive numerator is live. The oracle solves the concerted network of n protomers,
    # each moving through its conformation's forms, all flipping together at free E.
    mechs = (
        @allosteric_mechanism(begin
            substrates: A, B
            products: P
            catalytic_multiplicity: 1
            steps: begin
                E + A <--> E(A)        :: EqualAI
                E(A) + B ⇌ E(A, B)     :: EqualAI
                E(A, B) <--> E(P)      :: NonequalAI
                E + P ⇌ E(P)           :: EqualAI
            end
        end),
        @allosteric_mechanism(begin
            substrates: A, B
            products: P
            catalytic_multiplicity: 2
            steps: begin
                E + A <--> E(A)        :: EqualAI
                E(A) + B ⇌ E(A, B)     :: EqualAI
                E(A, B) <--> E(P)      :: NonequalAI
                E + P ⇌ E(P)           :: EqualAI
            end
        end),
        @allosteric_mechanism(begin
            substrates: A, B
            products: P
            catalytic_multiplicity: 3
            steps: begin
                E + A <--> E(A)        :: EqualAI
                E(A) + B ⇌ E(A, B)     :: EqualAI
                E(A, B) <--> E(P)      :: NonequalAI
                E + P ⇌ E(P)           :: EqualAI
            end
        end))
    for (n, em) in enumerate(mechs)
        am = _testhelper_lift(em)
        s = EnzymeRates.Derivation._solve(am)
        A = EnzymeRates.Derivation._law(am, s)
        @test A.n == n
        rng = Random.MersenneTwister(30 + n)
        for _ in 1:2
            point, concs, params = _testhelper_law_inputs(am, s.fitted, rng)
            oracle = setprecision(() -> _testhelper_protomer_rate(am, point, concs; n),
                                  BigFloat, 256)
            law = _testhelper_evaluate_law(A, merge(params, concs))
            @test abs(law - oracle) <= big(10)^-20 * abs(oracle)
        end
    end
end

@testset "an :OnlyI law is the exact limit of its vanishing constants" begin
    # The thirteen mechanisms whose :OnlyI groups bind modifiers alone, read without the
    # guard: each law equals the rapid-equilibrium limit of the oracle's own deletion
    # network (`_testhelper_limit_law`) at 3 points, and the active conformation's
    # denominator holds no constant of an :OnlyI group, since `conformation_steps` drops
    # them.
    for (i, (label, am)) in enumerate(pairs(_testhelper_onlyi_fixtures()))
        @testset "$label" begin
            sub = _testhelper_single_subunit(am)
            s = EnzymeRates.Derivation._solve(sub)
            A = EnzymeRates.Derivation._law(sub, s)
            rng = Random.MersenneTwister(400 + i)
            for _ in 1:3
                point, concs, params = _testhelper_law_inputs(sub, s.fitted, rng)
                exact = setprecision(BigFloat, 256) do
                    _testhelper_limit_law(sub, point, concs, Set{Symbol}(); rng)
                end
                law = _testhelper_evaluate_law(A, merge(params, concs))
                @test abs(law - exact) <= big(10)^-20 * abs(exact)
            end
            cm = EnzymeRates.Mechanisms.catalytic(sub)
            onlyi =
                Set(x for (k, ns) in enumerate(EnzymeRates.Derivation._edge_names(sub, :I))
                    if sub.mwc_tags[cm.steps[k].group] === :OnlyI
                    for x in ns)
            held = Set(A.vars.names[v] for mo in keys(A.active.D.terms) for (v, _) in mo)
            @test !isempty(onlyi)
            @test isempty(intersect(onlyi, held))
        end
    end
end

@testset "a flip-only inactive conformation is Q = 0, D = 1" begin
    # Both bindings :OnlyA: the inactive conformation's free E binds nothing, so that
    # conformation is free E alone. It turns nothing over and weighs 1, so the law is
    # H·QA/(DA + L), and it adds nothing to the term count.
    inert = @allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)      :: OnlyA
            E(S) <--> E(P)    :: EqualAI
            E + P ⇌ E(P)      :: OnlyA
        end
    end
    am = _testhelper_lift(inert)
    @test !any(EnzymeRates.Derivation.conformation_steps(am, :inactive))
    s = EnzymeRates.Derivation._solve(am)
    A = EnzymeRates.Derivation._law(am, s)
    @test EnzymeRates.Derivation._flip_only(A.inactive) &&
          !EnzymeRates.Derivation._flip_only(A.active)
    @test iszero(A.inactive.Q) && A.inactive.D == one(EnzymeRates.Derivation._Poly)
    @test A.inactive.trees == [one(EnzymeRates.Derivation._Poly)] &&
          isempty(A.inactive.shift)
    @test A.fA == one(EnzymeRates.Derivation._Poly) &&
          A.fI == one(EnzymeRates.Derivation._Poly)
    @test EnzymeRates.Derivation.rate_equation_complexity(am) == 1
    rng = Random.MersenneTwister(77)
    for _ in 1:3
        point, concs, params = _testhelper_law_inputs(am, s.fitted, rng)
        oracle = _testhelper_steady_state_rate(am, point, concs)
        law = _testhelper_evaluate_law(A, merge(params, concs))
        @test abs(law - oracle) <= big(10)^-20 * abs(oracle)
    end
end

@testset "_cross_weights divides out the common integer and monomial content" begin
    vars = EnzymeRates.Derivation._Vars()
    B = EnzymeRates.Derivation._var!(vars, :B; conc = true)
    k = EnzymeRates.Derivation._var!(vars, :k_ES_to_EP)
    q = EnzymeRates.Derivation._var!(vars, :K_EQ_to_E_Q)
    t(c, ps...) = EnzymeRates.Derivation._term(EnzymeRates.Derivation._monomial(ps...), c)
    # Two monomials, 2·k·B and 4·B², share the content 2·B.
    @test EnzymeRates.Derivation._cross_weights(t(2, k => 1, B => 1), t(4, B => 2)) ==
          (t(1, k => 1), t(2, B => 1))
    # Equal weights both become 1, whatever their polynomial part.
    w = t(1, k => 1) + t(1, B => 1)
    @test EnzymeRates.Derivation._cross_weights(w, w) ==
          (one(EnzymeRates.Derivation._Poly), one(EnzymeRates.Derivation._Poly))
    # A two-term weight keeps its primitive part; only the shared B comes out.
    @test EnzymeRates.Derivation._cross_weights(t(1, q => 1, B => 1) * w, t(1, B => 1)) ==
          (t(1, q => 1) * w, one(EnzymeRates.Derivation._Poly))
end

@testset "the cap counts both conformations" begin
    # A random-order bi-bi with the abortive complexes E(A, Q) and E(B, P), every step
    # steady state: 9 forms and 480 spanning trees, 4,320 products. With its chemistry
    # :NonequalAI at one catalytic subunit both conformations keep every step, so the law
    # counts 4,320 + 4,320 = 8,640: over the cap of 5,000, though each conformation is
    # under it. The engine refuses it before any polynomial product, with
    # RateEquationTooLargeError's message.
    em = @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)          :: EqualAI
            E + B <--> E(B)          :: EqualAI
            E(A) + B <--> E(A, B)    :: EqualAI
            E(B) + A <--> E(A, B)    :: EqualAI
            E(A, B) <--> E(P, Q)     :: NonequalAI
            E(P, Q) <--> E(P) + Q    :: EqualAI
            E(P, Q) <--> E(Q) + P    :: EqualAI
            E(P) <--> E + P          :: EqualAI
            E(Q) <--> E + Q          :: EqualAI
            E(A) + Q <--> E(A, Q)    :: EqualAI
            E(Q) + A <--> E(A, Q)    :: EqualAI
            E(B) + P <--> E(B, P)    :: EqualAI
            E(P) + B <--> E(B, P)    :: EqualAI
        end
    end
    am = _testhelper_lift(em)
    cm = am.catalytic
    amask, imask = EnzymeRates.Derivation.conformation_steps(am, :active),
        EnzymeRates.Derivation.conformation_steps(am, :inactive)
    @test _testhelper_laplacian_tree_count(cm, amask) == 4320
    @test _testhelper_laplacian_tree_count(cm, imask) == 4320
    @test EnzymeRates.Derivation.rate_equation_complexity(cm, amask) <=
          EnzymeRates.Derivation.MAX_RATE_EQUATION_TERMS
    @test EnzymeRates.Derivation.rate_equation_complexity(am) == 8640
    @test_throws(EnzymeRates.Derivation.RateEquationTooLargeError,
                 EnzymeRates.Derivation._law(am))
    @test_throws "has more than 5000 polynomial terms" EnzymeRates.Derivation._law(am)
    @test_throws EnzymeRates.Derivation.RateEquationTooLargeError rate_equation_string(em)
end
