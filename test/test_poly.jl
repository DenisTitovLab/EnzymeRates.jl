# ABOUTME: Tests of the exact sparse polynomials (src/poly.jl): monomials, checked
# ABOUTME: arithmetic, division by the Haldane factor, substitution, the shifts and content.

@testset "monomials multiply, raise and read exponents exactly" begin
    P = EnzymeRates
    @test P._monomial(2 => 1, 1 => 1//2, 2 => 2) == [1 => 1//2, 2 => 3//1]
    @test P._monomial(1 => 1, 1 => -1) == P._Mono()
    @test P._monomial(1 => 0, 3 => 2) == [3 => 2//1]
    @test P._mmul(P._monomial(1 => 1//2), P._monomial(1 => -1//2, 3 => 2)) == [3 => 2//1]
    @test P._mpow(P._monomial(1 => 1//2, 2 => -1), 3) == [1 => 3//2, 2 => -3//1]
    @test P._mpow(P._monomial(1 => 2), 1//2) == [1 => 1//1]
    @test P._mpow(P._monomial(1 => 2), 0) == P._Mono()
    @test P._mexp(P._monomial(1 => 1//2), 1) == 1//2
    @test P._mexp(P._monomial(1 => 1//2), 2) == 0
end

@testset "polynomials add, subtract and multiply exactly" begin
    P = EnzymeRates
    x, y = P._term(P._monomial(1 => 1)), P._term(P._monomial(2 => 1))
    p = x + y
    @test (p * p).terms == Dict(P._monomial(1 => 2) => 1, P._monomial(1 => 1, 2 => 1) => 2,
                                P._monomial(2 => 2) => 1)
    @test P._nterms(p * p) == 3
    @test (x - y) * (x + y) == x * x - y * y
    @test p - y == x && -x == P._scale(x, -1)
    @test iszero(p - p) && iszero(zero(P._Poly)) && !iszero(one(P._Poly))
    @test p * P._monomial(1 => -1) == one(P._Poly) + P._term(P._monomial(1 => -1, 2 => 1))
    @test P._scale(p, 0) == zero(P._Poly) && P._term(P._Mono(), 0) == zero(P._Poly)
    @test P._addto!(P._term(P._monomial(1 => 1), 2), P._monomial(1 => 1), -2) ==
          zero(P._Poly)
    # Equal by terms, hashed by terms, and distinct objects.
    @test x + y == y + x && hash(x + y) == hash(y + x) && x + y !== y + x
    @test ismutabletype(P._Poly) && isconst(P._Poly, :terms)
    @test fieldtypes(P._Poly) == (Dict{P._Mono, Int},)
end

@testset "an Int64 overflow throws OverflowError" begin
    P = EnzymeRates
    m = P._monomial(1 => 1)
    @test_throws OverflowError P._term(m, typemax(Int)) + P._term(m, 1)
    @test_throws OverflowError P._term(m, typemax(Int) ÷ 2 + 1) * P._term(P._Mono(), 2)
    @test_throws OverflowError P._scale(P._term(m, typemax(Int)), 2)
end

@testset "_divexact divides out the Haldane factor" begin
    P = EnzymeRates
    # Variables: 1 = A, 2 = P, 3 = Keq.
    mono(ps...) = P._term(P._monomial(ps...))
    H = mono(1 => 1) - mono(2 => 1, 3 => -1)
    q = mono(1 => 1, 3 => 1//2) + P._scale(mono(2 => 1, 3 => -1//2), 2) + one(P._Poly)
    @test P._divexact(H * q, H) == q
    # A substrate with two copies: H = A² − P/Keq.
    H2 = mono(1 => 2) - mono(2 => 1, 3 => -1)
    q2 = mono(1 => 1) + mono(2 => 1) + mono(1 => 1, 2 => 1, 3 => 1//2)
    @test P._divexact(H2 * q2, H2) == q2
    @test iszero(P._divexact(zero(P._Poly), H))
    @test_throws ArgumentError P._divexact(H * q + one(P._Poly), H)
    @test_throws ArgumentError P._divexact(q, q)
end

@testset "_subst replaces variables by power products" begin
    P = EnzymeRates
    mono(ps...) = P._term(P._monomial(ps...))
    p = mono(1 => 1, 2 => 2) + P._scale(mono(3 => 1), 3)
    f(v) = v == 2 ? P._monomial(4 => 1//2, 5 => -1) : nothing
    @test P._subst(f, p) == mono(1 => 1, 4 => 1, 5 => -2) + P._scale(mono(3 => 1), 3)
    # Terms that substitution makes alike merge, and cancel to zero.
    g(v) = v == 2 ? P._monomial(1 => 1) : nothing
    @test P._subst(g, mono(1 => 1) + mono(2 => 1)) == P._scale(mono(1 => 1), 2)
    @test iszero(P._subst(g, mono(1 => 1) - mono(2 => 1)))
end

@testset "the shifts bring polynomials to lowest terms in the concentrations" begin
    P = EnzymeRates
    t = P._Vars()
    a, b, k = P._var!(t, :A; conc = true), P._var!(t, :B; conc = true), P._var!(t, :k)
    @test (a, b, k) == (1, 2, 3) && P._var!(t, :A; conc = true) == 1
    @test t.names == [:A, :B, :k] && t.conc == [true, true, false]
    @test_throws ErrorException P._var!(t, :k; conc = true)
    @test t == P._Vars([:A, :B, :k], BitVector([1, 1, 0]), Dict(:A => 1, :B => 2, :k => 3))
    mono(ps...) = P._term(P._monomial(ps...))
    Q = mono(a => 1, b => 1, k => 1) + mono(a => 2, k => -1)
    D = mono(a => 1) + mono(a => 1, b => 1)
    # A is in every term, lowest 1; B is missing from some, lowest 0; k reaches −1.
    @test P._lowest_shift((Q, D)) == P._monomial(a => -1, k => 1)
    @test P._conc_shift((Q, D), t) == P._monomial(a => -1)
    D2 = mono(a => 1, k => 1) + mono(a => 1, b => 1) + P._scale(mono(k => 1//2), 2) +
         mono(b => 1, k => 1)
    @test P._split_conc(D2, t) == Dict(P._monomial(a => 1) => mono(k => 1),
                                      P._monomial(a => 1, b => 1) => one(P._Poly),
                                      P._Mono() => P._scale(mono(k => 1//2), 2),
                                      P._monomial(b => 1) => mono(k => 1))
end

@testset "_content splits off the integer and monomial content" begin
    P = EnzymeRates
    mono(ps...) = P._term(P._monomial(ps...))
    p = P._scale(mono(1 => 2, 3 => 1), 6) + P._scale(mono(1 => 1, 3 => 1//2), 4)
    @test P._content(p) ==
          (2, P._monomial(1 => 1, 3 => 1//2),
           P._scale(mono(1 => 1, 3 => 1//2), 3) + P._scale(one(P._Poly), 2))
    @test P._content(P._scale(mono(1 => 1), -2) + P._scale(one(P._Poly), -4)) ==
          (2, P._Mono(), -mono(1 => 1) - P._scale(one(P._Poly), 2))
    @test P._content(P._scale(mono(1 => -1, 2 => 1), 5)) ==
          (5, P._monomial(1 => -1, 2 => 1), one(P._Poly))
end
