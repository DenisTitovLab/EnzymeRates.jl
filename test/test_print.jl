# ABOUTME: Tests of the unit rule (src/normalize.jl) and the printed laws (src/print.jl):
# ABOUTME: each rule's normalizer, both printed forms against the law, names, MWC lines.

# The unit rule's normalizer of the plain mechanism `m`'s law, the law and the concrete
# mechanism.
function _testhelper_normalizer(m)
    mech = _testhelper_lift(m)
    L = EnzymeRates._law(mech)
    EnzymeRates._normalizer(L.D, L.vars, EnzymeRates._substrate_vars(mech, L.vars)), L, mech
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
    nz, L, mech = _testhelper_normalizer(mm)
    @test nz.rule === :constant
    @test nz.lambda == [EnzymeRates._Mono() => 1 // 1]
    @test isempty(nz.shift)
    _testhelper_check_units(mech, L)
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
    nz, L, mech = _testhelper_normalizer(pp)
    m(ps...) = _testhelper_conc_mono(L, ps...)
    @test nz.rule === :graded
    @test nz.lambda == [m(:A => 1) => 1 // 1, m(:B => 1) => 1 // 1,
                        m(:A => 1, :B => 1) => -1 // 1]
    @test isempty(nz.shift)
    _testhelper_check_units(mech, L)
end

@testset "Segel Bi Uni Uni Uni PP Ter Bi takes the graded basis with C" begin
    # The degree-1 term C enters the basis before Cleland's monomials A·B, A·C, B·C, so
    # the weights are C:2, A·B:1, A·C:−1, B·C:−1. The unit rule keeps this graded basis:
    # N is not Cleland's N = c(B·C)·c(A·C)·c(A·B)/c(A·B·C)².
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
    nz, L, mech = _testhelper_normalizer(terbi)
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
    _testhelper_check_units(mech, L)
end

@testset "the ping-pong power case puts A² before A·B" begin
    # Census ping-pong mechanism 6e46f55b0b68e5d9. A² and A·B are both degree 2 with no
    # non-substrate; the earlier name's larger exponent comes first, so A² joins A, B, P
    # and Q in the basis, and the weights are A:2, A²:−1.
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
    nz, L, mech = _testhelper_normalizer(power)
    m(ps...) = _testhelper_conc_mono(L, ps...)
    @test nz.rule === :graded
    @test nz.lambda == [m(:A => 1) => 2 // 1, m(:A => 2) => -1 // 1]
    _testhelper_check_units(mech, L)
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
    nz, L, mech = _testhelper_normalizer(f1)
    a = _testhelper_conc_mono(L, :A => 1)
    @test nz.rule === :ratio
    @test nz.lambda == [a => 1 // 1] && nz.shift == a
    _testhelper_check_units(mech, L)
end

@testset "F2 takes fractional weights" begin
    # F1 with F binding two A at rapid equilibrium: the terms A, P, A³ and A²·P. The basis
    # A, P, A³ solves the units equation with A:3/2, A³:−1/2.
    f2 = @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A <--> F + P
            F + A <--> E + P
            F + A + A ⇌ F(A, A)
        end
    end
    nz, L, mech = _testhelper_normalizer(f2)
    m(ps...) = _testhelper_conc_mono(L, ps...)
    @test nz.rule === :graded
    @test nz.lambda == [m(:A => 1) => 3 // 2, m(:A => 3) => -1 // 2]
    _testhelper_check_units(mech, L)
end

@testset "F3, a rapid-equilibrium Theorell–Chance step, divides by B" begin
    # Every term holds B or P (the rapid-equilibrium step E(A) + B ⇌ E(Q) + P ties them),
    # and the units equation has no solution on the basis B, P, A·B, B·Q, so N is the
    # first term, B with its coefficient.
    f3 = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B ⇌ E(Q) + P
            E(Q) <--> E + Q
        end
    end
    nz, L, mech = _testhelper_normalizer(f3)
    b = _testhelper_conc_mono(L, :B => 1)
    @test nz.rule === :ratio
    @test nz.lambda == [b => 1 // 1] && nz.shift == b
    _testhelper_check_units(mech, L)
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
    nz, L, mech = _testhelper_normalizer(f4)
    m(ps...) = _testhelper_conc_mono(L, ps...)
    @test nz.rule === :graded
    @test nz.lambda == [m(:A => 1) => 2 // 1, m(:A => 2) => -1 // 1]
    _testhelper_check_units(mech, L)
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
        subs = EnzymeRates._substrate_vars(am, A.vars)
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
        # The normalized polynomials in the MWC combination give the law. That holds for
        # any nonzero N_A and N_I, so it pins the invariance; the rule and weight pins
        # above test the normalizer itself.
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
    # Units on the 42 specs, each at one catalytic subunit (spec section 14).
    for spec in MECHANISM_TEST_SPECS
        @testset "$(spec.name)" begin
            am = _testhelper_lift(_testhelper_single_subunit(spec.mechanism))
            _testhelper_check_units(am, EnzymeRates._law(am, EnzymeRates._solve(am)))
        end
    end
end

@testset "the printed law and its coefficient form equal the law on every spec" begin
    # Both forms, read back in 256-bit BigFloat, equal the derived law to 1e-20 at a
    # consistent point of every spec, the concrete and the compiled mechanism print one
    # text, and neither form holds a printing artifact.
    for spec in MECHANISM_TEST_SPECS
        @testset "$(spec.name)" begin
            m = spec.mechanism
            full = rate_equation_string(m)
            @test rate_equation_string(_testhelper_lift(m)) == full
            @test isempty(_testhelper_print_artifacts(full))
            @test isempty(_testhelper_print_artifacts(
                rate_equation_string(m; coefficients = true)))
            @test _testhelper_printed_law_error(m; points = 1) <= big(10)^-20
        end
    end
end

@testset "coefficient names read powers and ratios, and N prints the normalizer once" begin
    # F1 reads A and P only through P/A, so the ratio form divides by A: the coefficients
    # of P/A take ⁻¹, and the polynomial factor of the divisor prints as N.
    f1 = @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A <--> F + P
            F + A <--> E + P
        end
    end
    c1 = split(rate_equation_string(f1; coefficients = true), '\n')
    @test c1[1] == "v = E_total * (n_0 + n_P_A⁻¹ * P / A) * (A - P / Keq) / " *
                   "(1 + d_P_A⁻¹ * P / A)"
    @test c1[2] == "# The coefficients are functions of the fitted constants, not " *
                   "parameters:"
    @test c1[3:end] == [
        "N = k_E_A_to_F_P + k_F_A_to_E_P",
        "n_0 = 2 * k_E_A_to_F_P * k_F_A_to_E_P / N",
        "n_P_A⁻¹ = 2 * k_E_A_to_F_P * k_F_A_to_E_P / (Keq * N)",
        "d_P_A⁻¹ = (k_E_A_to_F_P * k_F_A_to_E_P + k_F_P_to_E_A^2 * Keq^2) / " *
        "(k_F_P_to_E_A * Keq^2 * N)"]
    # F3's rapid-equilibrium Theorell–Chance step makes every denominator term carry B or
    # P: the ratio form divides by B, whose coefficient is one constant, so no N prints.
    f3 = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B ⇌ E(Q) + P
            E(Q) <--> E + Q
        end
    end
    c3 = split(rate_equation_string(f3; coefficients = true), '\n')
    @test c3[1] == "v = E_total * (n_B⁻¹ / B) * (A * B - P * Q / Keq) / (1 + " *
                   "d_P_B⁻¹ * P / B + d_A * A + d_A_P_B⁻¹ * A * P / B + d_Q * Q + " *
                   "d_P_Q_B⁻¹ * P * Q / B)"
    @test "n_B⁻¹ = k_E_A_to_EA" in c3
    @test !any(startswith("N = "), c3)
    # F2's weights are fractional: N is a power 3/2, and its coefficients hold roots.
    f2 = @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A <--> F + P
            F + A <--> E + P
            F + A + A ⇌ F(A, A)
        end
    end
    @test occursin("\nN = (k_E_A_to_F_P + k_F_A_to_E_P)^(3/2)\n", rate_equation_string(f2))
    @test "n_A = 2 * sqrt(K_F_A_A_to_FAA) * k_E_A_to_F_P^(3/2) * k_F_A_to_E_P / N" in
          split(rate_equation_string(f2; coefficients = true), '\n')
    # F2 with its two steps in one group: the weight 3/2 falls on the coefficient of A,
    # whose integer content is 2. An integer to an integer power divides the coefficients;
    # to the power 3/2 it prints apart, as N.
    f2g = @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            (E + A <--> F + P, F + A <--> E + P)
            F + A + A ⇌ F(A, A)
        end
    end
    @test occursin("\nN = (2)^(3/2)\n", rate_equation_string(f2g))
    # The random-order bi-bi's denominator holds A^2 * B, named by its superscript.
    ro = only(s for s in MECHANISM_TEST_SPECS if s.name == "Random-order Bi-Bi").mechanism
    @test occursin(" + d_A²_B * A^2 * B + ", rate_equation_string(ro; coefficients = true))
    # The unit-rule testsets above pin the weights of F4, Segel's ter bi and the ping-pong
    # power case; N prints the normalizer's polynomial factor alone, its integer
    # and monomial content moved into the coefficients. F4's is a sum squared. Segel's
    # graded N is the product of two sums, where Cleland's would divide it by a third sum
    # squared. The power case's is one sum squared, Cleland's N squared.
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
    @test occursin("\nN = (k_E_A_to_F_P + k_F_A_to_E_P)^2\n", rate_equation_string(f4))
    terbi = only(s for s in MECHANISM_TEST_SPECS
                 if s.name == "Segel Bi Uni Uni Uni PP Ter Bi").mechanism
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
    n_line(m) = only(filter(startswith("N = "), split(rate_equation_string(m), '\n')))
    @test occursin(r"^N = \([^()]+\) \* \([^()]+\)$", n_line(terbi))
    @test occursin(r"^N = \([^()]+\)\^2$", n_line(power))
    for m in (f1, f2, f2g, f3, f4, terbi, power)
        @test _testhelper_printed_law_error(m) <= big(10)^-20
        @test isempty(_testhelper_print_artifacts(rate_equation_string(m))) &&
              isempty(_testhelper_print_artifacts(
                  rate_equation_string(m; coefficients = true)))
    end
end

@testset "an :OnlyI inhibitor prints a law equal to its limit" begin
    # I binds only the inactive conformation, at rapid equilibrium and at
    # steady state. Both printed forms equal the rapid-equilibrium limit of the oracle's
    # own network (`_testhelper_limit_law`) to 1e-20 at two points; the active
    # conformation's denominator holds no constant of I's binding, the inactive one's does.
    competitive_RE = @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)            :: EqualAI
            E(S) <--> E + P         :: EqualAI
            E + I ⇌ E(I)            :: OnlyI
        end
    end
    competitive_SS = @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)            :: EqualAI
            E(S) <--> E + P         :: NonequalAI
            E + I <--> E(I)         :: OnlyI
        end
    end
    cases = ((competitive_RE, "D_I = 1 + S / K_ES_to_E_S + I / K_I_EI_to_E_I",
              "dI_I = 1 / K_I_EI_to_E_I"),
             (competitive_SS,
              "D_I = 1 + S / K_ES_to_E_S + I * k_I_E_I_to_EI / k_I_EI_to_E_I",
              "dI_I = k_I_E_I_to_EI / k_I_EI_to_E_I"))
    for (i, (m, d_i, d_i_coefficient)) in enumerate(cases)
        full = rate_equation_string(m)
        coef = rate_equation_string(m; coefficients = true)
        @test "D_A = 1 + S / K_ES_to_E_S" in split(full, '\n')
        @test d_i in split(full, '\n')
        @test last(split(full, '\n')) ==
              "v = E_total * (Q_A + L * Q_I) * (S - P / Keq) / (D_A + L * D_I)"
        @test "D_A = 1 + dA_S * S" in split(coef, '\n')
        @test d_i_coefficient in split(coef, '\n')
        am = _testhelper_lift(m)
        cm = EnzymeRates.catalytic(am)
        onlyi = [String(x) for (k, ns) in enumerate(EnzymeRates._edge_names(am, :I))
                 if EnzymeRates.cat_allo_state(am, cm.group[k]) === :OnlyI for x in ns]
        d_a = only(l for l in split(full, '\n') if startswith(l, "D_A = "))
        @test !isempty(onlyi)
        @test !any(x -> occursin(x, d_a), onlyi)
        @test any(x -> occursin(x, d_i), onlyi)
        s = EnzymeRates._solve(am)
        rng = Random.MersenneTwister(500 + i)
        for _ in 1:2
            point, concs, params = _testhelper_law_inputs(am, s.fitted, rng)
            exact = setprecision(BigFloat, 256) do
                _testhelper_limit_law(am, point, concs, Set{Symbol}(); rng)
            end
            p = merge(params, Dict(:E_total => big(1)))
            for text in (full, coef)
                @test abs(_testhelper_eval_rate_string(text, p, concs) - exact) <=
                      big(10)^-20 * abs(exact)
            end
        end
        @test _testhelper_printed_law_error(m) <= big(10)^-20
    end
end

@testset "the MWC lines name each conformation's polynomials" begin
    # A two-site homodimer prints Q, D and W for each conformation and no F: both
    # free-enzyme weights are 1. test_rate_eq_derivation.jl pins its text.
    two_site = @allosteric_mechanism begin
        substrates: S
        products:   P
        allosteric_regulators: R
        catalytic_multiplicity: 2
        steps: begin
            E + P ⇌ E(P)     :: NonequalAI
            E + S ⇌ E(S)     :: NonequalAI
            E(S) <--> E(P)   :: NonequalAI
        end
        regulatory_site(2): R::NonequalAI
    end
    # A ping-pong whose second half runs in the active conformation only: neither
    # conformation's denominator has a constant term, so each prints its own normalizer,
    # and free E weighs F6P in one and F16BP in the other, so F_A and F_I print.
    half2 = @allosteric_mechanism begin
        substrates: ATP, F6P
        products: ADP, F16BP
        catalytic_multiplicity: 1
        steps: begin
            E + ATP ⇌ E(ATP)                                                      :: EqualAI
            E(ATP) <--> E(F16BP; residual = ATP - F16BP)                          :: EqualAI
            E(; residual = ATP - F16BP) + F16BP ⇌ E(F16BP; residual = ATP - F16BP):: EqualAI
            E(; residual = ATP - F16BP) + F6P ⇌ E(F6P; residual = ATP - F16BP)    :: OnlyA
            E(F6P; residual = ATP - F16BP) ⇌ E(ADP)                               :: OnlyA
            E + ADP ⇌ E(ADP)                                                      :: EqualAI
        end
    end
    # Both bindings :OnlyA: the inactive conformation is free E alone, so the
    # denominator's inactive term is L with no factor.
    inert = @allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)      :: OnlyA
            E(S) <--> E(P)    :: EqualAI
            E + P ⇌ E(P)      :: OnlyA
        end
    end
    heads(t) = [l[1:3] for l in split(t, '\n') if occursin(r"^[QDFWN]_[AI] = ", l)]
    @test heads(rate_equation_string(two_site)) ==
          ["Q_A", "D_A", "Q_I", "D_I", "W_A", "W_I"]
    t2 = rate_equation_string(half2)
    @test heads(t2) == ["N_A", "Q_A", "D_A", "N_I", "D_I", "F_A", "F_I"]
    @test "F_A = F6P / (K_A_EF6P_res_plusATP_minusF16BP_to_E_res_plusATP_minusF16BP_F6P " *
          "* N_A)" in split(t2, '\n')
    @test last(split(t2, '\n')) == "v = E_total * (F_I * Q_A) * (ATP * F6P - ADP * " *
                                    "F16BP / Keq) / (F_I * D_A + L * F_A * D_I)"
    ti = rate_equation_string(inert)
    @test heads(ti) == ["Q_A", "D_A"]
    @test last(split(ti, '\n')) == "v = E_total * (Q_A * D_A) * (S - P / Keq) / (D_A^2 + L)"
    for m in (two_site, half2, inert)
        @test _testhelper_printed_law_error(m) <= big(10)^-20
        @test isempty(_testhelper_print_artifacts(rate_equation_string(m))) &&
              isempty(_testhelper_print_artifacts(
                  rate_equation_string(m; coefficients = true)))
    end
end

@testset "the printed-law reader runs the full law only as written" begin
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
        end
    end
    s = rate_equation_string(m)
    p = _testhelper_random_independent_params(m; rng = Random.MersenneTwister(1))
    c = (S = 0.4, P = 0.9)
    @test isapprox(_testhelper_eval_rate_string(s, p, c), rate_equation(m, c, p);
                   rtol = 1e-12)
    # A name the params or concs line leaves out, and a line read before the lines that
    # bind its names, fail as they do in Julia.
    for t in (replace(s, "K_ES_to_E_S, " => ""; count = 1),
              replace(s, "(; S, P) = concs" => "(; S) = concs"),
              join(reverse(split(s, '\n')), '\n'))
        @test_throws "before it binds it" _testhelper_eval_rate_string(t, p, c)
    end
end
