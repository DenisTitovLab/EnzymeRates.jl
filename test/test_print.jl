# ABOUTME: Tests of what a printed rate law rests on: the unit rule's normalizer of each
# ABOUTME: conformation's denominator, and the units of every derived law from the graph.

# The unit rule's normalizer of the plain mechanism `m`'s law, and the law.
function _testhelper_normalizer(m)
    mech = _testhelper_lift(m)
    L = EnzymeRates._law(mech)
    subs = Int[L.vars.index[x] for x in unique(EnzymeRates.substrates(mech))]
    EnzymeRates._normalizer(L.D, L.vars, subs), L
end

# The concentration monomial of law `L` with the exponents `x => e`.
_testhelper_conc_mono(L, ps::Pair{Symbol}...) =
    EnzymeRates._monomial((L.vars.index[x] => e for (x, e) in ps)...)

# ── The unit rule (spec section 9) ───────────────────────────────────────────

@testset "Michaelis–Menten takes its constant term" begin
    # Rule 1: the denominator has a term free of concentrations, the free enzyme's, and
    # N is its coefficient, so every term reads as a form's abundance relative to free E.
    mm = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> E(S)
            E(S) <--> E + P
        end
    end
    nz, L = _testhelper_normalizer(mm)
    @test nz.rule === :constant
    @test nz.lambda == [EnzymeRates._Mono() => 1 // 1]
    @test isempty(nz.shift)
    _testhelper_check_units(_testhelper_lift(mm), L)
end

@testset "Segel Ping Pong Bi Bi takes the graded basis A, B, A·B" begin
    # Rule 2: no constant term. A, B, P and Q (degree 1), then A·B, are affinely
    # independent; the weights A:1, B:1, A·B:−1 solve Σλ = 1 and Σλ·a = 0, and P and Q
    # weigh 0. Here the graded basis is Cleland's N: ΠS/X over each substrate X, over ΠS.
    pp = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) <--> F + P
            F + B <--> F(B)
            F(B) <--> E + Q
        end
    end
    nz, L = _testhelper_normalizer(pp)
    m(ps...) = _testhelper_conc_mono(L, ps...)
    @test nz.rule === :graded
    @test nz.lambda == [m(:A => 1) => 1 // 1, m(:B => 1) => 1 // 1,
                        m(:A => 1, :B => 1) => -1 // 1]
    @test isempty(nz.shift)
    _testhelper_check_units(_testhelper_lift(pp), L)
end

@testset "Segel Bi Uni Uni Uni PP Ter Bi takes the graded basis with C" begin
    # The degree-1 term C enters the basis before Cleland's monomials A·B, A·C, B·C, so
    # the weights are C:2, A·B:1, A·C:−1, B·C:−1. P5-10a keeps this graded basis: N is
    # not Cleland's N = c(B·C)·c(A·C)·c(A·B)/c(A·B·C)².
    terbi = @enzyme_mechanism begin
        substrates: A, B, C
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B <--> E(A, B)
            E(A, B) <--> F + P
            F + C <--> F(C)
            F(C) <--> E + Q
        end
    end
    nz, L = _testhelper_normalizer(terbi)
    m(ps...) = _testhelper_conc_mono(L, ps...)
    @test nz.rule === :graded
    @test nz.lambda == [m(:C => 1) => 2 // 1, m(:A => 1, :B => 1) => 1 // 1,
                        m(:A => 1, :C => 1) => -1 // 1, m(:B => 1, :C => 1) => -1 // 1]
    rng = Random.MersenneTwister(31)
    env = Dict{Symbol, BigFloat}(x => exp(2 * big(rand(rng)) - 1) for x in L.vars.names)
    Dc = EnzymeRates._split_conc(L.D, L.vars)
    c(ps...) = _testhelper_poly_value(Dc[m(ps...)], L.vars, env)
    cleland = c(:B => 1, :C => 1) * c(:A => 1, :C => 1) * c(:A => 1, :B => 1) /
              c(:A => 1, :B => 1, :C => 1)^2
    @test !isapprox(_testhelper_norm_value(nz, L.D, L.vars, env), cleland; rtol = 1e-6)
    _testhelper_check_units(_testhelper_lift(terbi), L)
end

@testset "the ping-pong power case puts A² before A·B" begin
    # Census ping-pong mechanism 6e46f55b0b68e5d9 (UNIT §2). A² and A·B are both degree
    # 2 with no non-substrate; the earlier name's larger exponent comes first (P5-11), so
    # A² joins A, B, P and Q in the basis, and the weights are A:2, A²:−1.
    power = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A <--> E(P; residual = A - P), E(Q) + A <--> E(A, Q))
            (E + Q <--> E(Q), E(P; residual = A - P) + Q <--> E(A, Q))
            (E(; residual = A - P) + B <--> E(Q),
             E(P; residual = A - P) + B <--> E(B, P; residual = A - P))
            (E(; residual = A - P) + P ⇌ E(P; residual = A - P),
             E(Q) + P ⇌ E(B, P; residual = A - P))
        end
    end
    nz, L = _testhelper_normalizer(power)
    m(ps...) = _testhelper_conc_mono(L, ps...)
    @test nz.rule === :graded
    @test nz.lambda == [m(:A => 1) => 2 // 1, m(:A => 2) => -1 // 1]
    _testhelper_check_units(_testhelper_lift(power), L)
end

@testset "F1 divides by its first term, A" begin
    # Rule 3: two half-reactions between the free forms E and F, each A → P. Every term
    # holds A or P, and no weights on A and P solve the units equation, so the law reads A
    # and P only through their ratio: N is the first term, its monomial A included.
    f1 = @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A <--> F + P
            F + A <--> E + P
        end
    end
    nz, L = _testhelper_normalizer(f1)
    a = _testhelper_conc_mono(L, :A => 1)
    @test nz.rule === :ratio
    @test nz.lambda == [a => 1 // 1] && nz.shift == a
    _testhelper_check_units(_testhelper_lift(f1), L)
end

@testset "F2 takes fractional weights" begin
    # F1 with F binding two A at rapid equilibrium: the terms A, P, A³ and A²·P. The basis
    # A, P, A³ solves the units equation with A:3/2, A³:−1/2 (P5-12).
    f2 = @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A <--> F + P
            F + A <--> E + P
            F + A + A ⇌ F(A, A)
        end
    end
    nz, L = _testhelper_normalizer(f2)
    m(ps...) = _testhelper_conc_mono(L, ps...)
    @test nz.rule === :graded
    @test nz.lambda == [m(:A => 1) => 3 // 2, m(:A => 3) => -1 // 2]
    _testhelper_check_units(_testhelper_lift(f2), L)
end

@testset "F3, a rapid-equilibrium Theorell–Chance step, divides by B" begin
    # Every term holds B or P (the rapid-equilibrium step E(A) + B ⇌ E(Q) + P ties them),
    # and the units equation has no solution on the basis B, P, A·B, B·Q, so N is the
    # first term, B with its coefficient (P5-13).
    f3 = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B ⇌ E(Q) + P
            E(Q) <--> E + Q
        end
    end
    nz, L = _testhelper_normalizer(f3)
    b = _testhelper_conc_mono(L, :B => 1)
    @test nz.rule === :ratio
    @test nz.lambda == [b => 1 // 1] && nz.shift == b
    _testhelper_check_units(_testhelper_lift(f3), L)
end

@testset "F4 takes integer weights A:2, A²:−1" begin
    # F1 with F binding A twice in two rapid-equilibrium steps: the terms A, P, A², A·P,
    # A³ and A²·P. The basis A, P, A² solves the units equation with A:2, A²:−1.
    f4 = @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A <--> F + P
            F + A <--> E + P
            F + A ⇌ F(A)
            F(A) + A ⇌ F(A, A)
        end
    end
    nz, L = _testhelper_normalizer(f4)
    m(ps...) = _testhelper_conc_mono(L, ps...)
    @test nz.rule === :graded
    @test nz.lambda == [m(:A => 1) => 2 // 1, m(:A => 2) => -1 // 1]
    _testhelper_check_units(_testhelper_lift(f4), L)
end

@testset "each MWC conformation takes its own rule" begin
    # Spec section 9: an MWC law normalizes each conformation's polynomials by that
    # conformation's rule, dividing its Q, D and free-E weight alike, so the MWC
    # combination does not change. A two-subunit ping-pong whose chemistry differs between
    # the conformations takes Segel's graded basis in both. F3 with A's binding and the
    # Theorell–Chance step :OnlyA takes the ratio form in the active conformation and the
    # constant term in the inactive one, which keeps free E and E(Q).
    pp2 = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        catalytic_multiplicity: 2
        steps: begin
            E + A <--> E(A)       :: EqualAI
            E(A) <--> F + P       :: NonequalAI
            F + B <--> F(B)       :: EqualAI
            F(B) <--> E + Q       :: EqualAI
        end
    end)
    tc = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)            :: OnlyA
            E(A) + B ⇌ E(Q) + P        :: OnlyA
            E(Q) <--> E + Q            :: EqualAI
        end
    end)
    for (i, am) in enumerate((pp2, tc))
        s = EnzymeRates._solve(am)
        A = EnzymeRates._law(am, s)
        subs = Int[A.vars.index[x] for x in unique(EnzymeRates.substrates(am))]
        nA = EnzymeRates._normalizer(A.active.D, A.active.vars, subs)
        nI = EnzymeRates._normalizer(A.inactive.D, A.inactive.vars, subs)
        m(ps...) = _testhelper_conc_mono(A.active, ps...)
        if i == 1
            segel = [m(:A => 1) => 1 // 1, m(:B => 1) => 1 // 1,
                     m(:A => 1, :B => 1) => -1 // 1]
            @test (nA.rule, nI.rule) == (:graded, :graded)
            @test nA.lambda == nI.lambda == segel
        else
            @test (nA.rule, nI.rule) == (:ratio, :constant)
            @test nA.lambda == [m(:B => 1) => 1 // 1] && nA.shift == m(:B => 1)
            @test nI.lambda == [EnzymeRates._Mono() => 1 // 1]
        end
        _testhelper_check_units(am, A)
        rng = Random.MersenneTwister(32 + i)
        for _ in 1:3
            _, concs, params = _testhelper_law_inputs(am, s.fitted, rng)
            env = merge(params, concs)
            value(p) = _testhelper_poly_value(p, A.vars, env)
            NA = _testhelper_norm_value(nA, A.active.D, A.active.vars, env)
            NI = _testhelper_norm_value(nI, A.inactive.D, A.inactive.vars, env)
            n, l = A.n, env[:L]
            qa, da, fa = value(A.active.Q) / NA, value(A.active.D) / NA, value(A.fA) / NA
            qi, di, fi = (value(A.inactive.Q) / NI, value(A.inactive.D) / NI,
                          value(A.fI) / NI)
            num = fi^n * qa * da^(n - 1) + l * fa^n * qi * di^(n - 1)
            normalized = value(A.H) * num / (fi^n * da^n + l * fa^n * di^n)
            law = _testhelper_evaluate_law(A, env)
            @test abs(normalized - law) <= big(10)^-20 * abs(law)
        end
    end
end

@testset "every spec's law has one unit per denominator and a rate per time" begin
    # Q1's units on the 42 specs, each at one catalytic subunit (spec section 14).
    for spec in MECHANISM_TEST_SPECS
        @testset "$(spec.name)" begin
            am = _testhelper_lift(_testhelper_single_subunit(spec.mechanism))
            _testhelper_check_units(am, EnzymeRates._law(am, EnzymeRates._solve(am)))
        end
    end
end
