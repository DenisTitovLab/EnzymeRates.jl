# ABOUTME: Tests of the sparse polynomials (src/Derivation/poly.jl): monomials, checked
# ABOUTME: arithmetic, division by the Haldane factor, substitution, the shifts and content.

@testset "monomials multiply, raise and read exponents exactly" begin
    P = EnzymeRates
    @test P.Derivation._monomial(2 => 1, 1 => 1//2, 2 => 2) == [1 => 1//2, 2 => 3//1]
    @test P.Derivation._monomial(1 => 1, 1 => -1) == P.Derivation._Mono()
    @test P.Derivation._monomial(1 => 0, 3 => 2) == [3 => 2//1]
    @test P.Derivation._mmul(P.Derivation._monomial(1 => 1//2),
                             P.Derivation._monomial(1 => -1//2, 3 => 2)) == [3 => 2//1]
    @test P.Derivation._mpow(P.Derivation._monomial(1 => 1//2, 2 => -1), 3) ==
          [1 => 3//2, 2 => -3//1]
    @test P.Derivation._mpow(P.Derivation._monomial(1 => 2), 1//2) == [1 => 1//1]
    @test P.Derivation._mpow(P.Derivation._monomial(1 => 2), 0) == P.Derivation._Mono()
    @test P.Derivation._mexp(P.Derivation._monomial(1 => 1//2), 1) == 1//2
    @test P.Derivation._mexp(P.Derivation._monomial(1 => 1//2), 2) == 0
end

@testset "polynomials add, subtract and multiply exactly" begin
    P = EnzymeRates
    x, y = P.Derivation._term(P.Derivation._monomial(1 => 1)),
        P.Derivation._term(P.Derivation._monomial(2 => 1))
    p = x + y
    @test (p * p).terms == Dict(P.Derivation._monomial(1 => 2) => 1,
                                P.Derivation._monomial(1 => 1, 2 => 1) => 2,
                                P.Derivation._monomial(2 => 2) => 1)
    @test P.Derivation._nterms(p * p) == 3
    @test (x - y) * (x + y) == x * x - y * y
    @test p - y == x && -x == P.Derivation._scale(x, -1)
    @test iszero(p - p) && iszero(zero(P.Derivation._Poly)) &&
          !iszero(one(P.Derivation._Poly))
    @test p * P.Derivation._monomial(1 => -1) == one(P.Derivation._Poly) +
          P.Derivation._term(P.Derivation._monomial(1 => -1, 2 => 1))
    @test P.Derivation._scale(p, 0) == zero(P.Derivation._Poly) &&
          P.Derivation._term(P.Derivation._Mono(), 0) == zero(P.Derivation._Poly)
    @test P.Derivation._addto!(P.Derivation._term(P.Derivation._monomial(1 => 1), 2),
                               P.Derivation._monomial(1 => 1), -2) ==
          zero(P.Derivation._Poly)
    # Equal by terms, hashed by terms, and distinct objects.
    @test x + y == y + x && hash(x + y) == hash(y + x) && x + y !== y + x
    @test ismutabletype(P.Derivation._Poly) && isconst(P.Derivation._Poly, :terms)
    @test fieldtypes(P.Derivation._Poly) == (Dict{P.Derivation._Mono, Int},)
end

@testset "an Int64 overflow throws OverflowError" begin
    P = EnzymeRates
    m = P.Derivation._monomial(1 => 1)
    @test_throws OverflowError P.Derivation._term(m, typemax(Int)) +
        P.Derivation._term(m, 1)
    @test_throws OverflowError P.Derivation._term(m, typemax(Int) ÷ 2 + 1) *
        P.Derivation._term(P.Derivation._Mono(), 2)
    @test_throws OverflowError P.Derivation._scale(P.Derivation._term(m, typemax(Int)), 2)
end

@testset "_divexact divides out the Haldane factor" begin
    P = EnzymeRates
    # Variables: 1 = A, 2 = P, 3 = Keq.
    mono(ps...) = P.Derivation._term(P.Derivation._monomial(ps...))
    H = mono(1 => 1) - mono(2 => 1, 3 => -1)
    q = mono(1 => 1, 3 => 1//2) + P.Derivation._scale(mono(2 => 1, 3 => -1//2), 2) +
        one(P.Derivation._Poly)
    @test P.Derivation._divexact(H * q, H) == q
    # A substrate with two copies: H = A² − P/Keq.
    H2 = mono(1 => 2) - mono(2 => 1, 3 => -1)
    q2 = mono(1 => 1) + mono(2 => 1) + mono(1 => 1, 2 => 1, 3 => 1//2)
    @test P.Derivation._divexact(H2 * q2, H2) == q2
    # A chain with a missing power: H·(A³ + P³/Keq³) has no A²·P²/Keq² term.
    qg = mono(1 => 3) + mono(2 => 3, 3 => -3)
    @test P.Derivation._divexact(H * qg, H) == qg
    @test iszero(P.Derivation._divexact(zero(P.Derivation._Poly), H))
    @test_throws ArgumentError P.Derivation._divexact(H * q + one(P.Derivation._Poly), H)
    @test_throws ArgumentError P.Derivation._divexact(q, q)
end

@testset "_subst replaces variables by power products" begin
    P = EnzymeRates
    mono(ps...) = P.Derivation._term(P.Derivation._monomial(ps...))
    p = mono(1 => 1, 2 => 2) + P.Derivation._scale(mono(3 => 1), 3)
    f(v) = v == 2 ? P.Derivation._monomial(4 => 1//2, 5 => -1) : nothing
    @test P.Derivation._subst(f, p) == mono(1 => 1, 4 => 1, 5 => -2) +
          P.Derivation._scale(mono(3 => 1), 3)
    # Terms that substitution makes alike merge, and cancel to zero.
    g(v) = v == 2 ? P.Derivation._monomial(1 => 1) : nothing
    @test P.Derivation._subst(g, mono(1 => 1) + mono(2 => 1)) ==
          P.Derivation._scale(mono(1 => 1), 2)
    @test iszero(P.Derivation._subst(g, mono(1 => 1) - mono(2 => 1)))
end

@testset "the shifts bring polynomials to lowest terms in the concentrations" begin
    P = EnzymeRates
    t = P.Derivation._Vars()
    a, b, k = P.Derivation._var!(t, :A; conc = true),
        P.Derivation._var!(t, :B; conc = true), P.Derivation._var!(t, :k)
    @test (a, b, k) == (1, 2, 3) && P.Derivation._var!(t, :A; conc = true) == 1
    @test t.names == [:A, :B, :k] && t.conc == [true, true, false]
    @test_throws ErrorException P.Derivation._var!(t, :k; conc = true)
    @test_throws ErrorException("internal error: A is a concentration in one use and a " *
                                "constant in another") P.Derivation._var!(t, :A)
    @test t == P.Derivation._Vars([:A, :B, :k], BitVector([1, 1, 0]),
                                  Dict(:A => 1, :B => 2, :k => 3))
    mono(ps...) = P.Derivation._term(P.Derivation._monomial(ps...))
    Q = mono(a => 1, b => 1, k => 1) + mono(a => 2, k => -1)
    D = mono(a => 1) + mono(a => 1, b => 1)
    # A is in every term, lowest 1; B is missing from some, lowest 0; k reaches −1.
    @test P.Derivation._lowest_shift((Q, D)) == P.Derivation._monomial(a => -1, k => 1)
    @test P.Derivation._conc_shift((Q, D), t) == P.Derivation._monomial(a => -1)
    D2 = mono(a => 1, k => 1) + mono(a => 1, b => 1) +
        P.Derivation._scale(mono(k => 1//2), 2) + mono(b => 1, k => 1)
    @test P.Derivation._split_conc(D2, t) ==
          Dict(P.Derivation._monomial(a => 1) => mono(k => 1),
               P.Derivation._monomial(a => 1, b => 1) => one(P.Derivation._Poly),
               P.Derivation._Mono() => P.Derivation._scale(mono(k => 1//2), 2),
               P.Derivation._monomial(b => 1) => mono(k => 1))
end

@testset "_content splits off the integer and monomial content" begin
    P = EnzymeRates
    mono(ps...) = P.Derivation._term(P.Derivation._monomial(ps...))
    p = P.Derivation._scale(mono(1 => 2, 3 => 1), 6) +
        P.Derivation._scale(mono(1 => 1, 3 => 1//2), 4)
    @test P.Derivation._content(p) == (2, P.Derivation._monomial(1 => 1, 3 => 1//2),
                                       P.Derivation._scale(mono(1 => 1, 3 => 1//2), 3) +
                                       P.Derivation._scale(one(P.Derivation._Poly), 2))
    @test P.Derivation._content(P.Derivation._scale(mono(1 => 1), -2) +
                                P.Derivation._scale(one(P.Derivation._Poly), -4)) ==
          (2, P.Derivation._Mono(),
           -mono(1 => 1) - P.Derivation._scale(one(P.Derivation._Poly), 2))
    @test P.Derivation._content(P.Derivation._scale(mono(1 => -1, 2 => 1), 5)) ==
          (5, P.Derivation._monomial(1 => -1, 2 => 1), one(P.Derivation._Poly))
end
