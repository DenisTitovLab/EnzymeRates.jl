# ABOUTME: Checks poly.jl on random polynomials: ring identities, evaluation with rational
# ABOUTME: exponents in BigFloat, substitution, exact binomial division, lowest terms, overflow.
include(joinpath(@__DIR__, "..", "src", "Slice.jl"))
using Random, Printf
const S = Slice
setprecision(BigFloat, 512)

rng = Xoshiro(7)
vars = S.Vars()
conc = [S.var!(vars, x; conc = true) for x in (:A, :B, :P, :Q)]
cons = [S.var!(vars, Symbol("k", i)) for i in 1:5]
keq = S.var!(vars, :Keq)
randexp(c) = c ? rand(rng, -1:2) // 1 : rand(rng, (-1//2, 1//2, 1//1, -1//1, 3//2, 0//1))
function randpoly(T; n = rand(rng, 1:6), laurent = false)
    p = S.Poly{T}()
    for _ in 1:n
        m = S.mono([v => (laurent ? randexp(true) : rand(rng, 0:2)) for v in conc]...,
                   [v => randexp(false) for v in rand(rng, cons, 2)]...)
        S.addto!(p, m, rand(rng, -5:5))
    end
    p
end
vals() = BigFloat[exp(BigFloat(randn(rng))) for _ in vars.names]
relerr(a, b) = abs(a - b) / max(abs(b), big(1e-300))

worst = Dict{String, BigFloat}()
note(k, e) = (worst[k] = max(get(worst, k, big(0)), e))
fails = String[]
for trial in 1:2000
    p, q, r = randpoly(Int), randpoly(Int), randpoly(Int)
    (p + q) * r == p * r + q * r || push!(fails, "distributivity $trial")
    p - p == zero(S.Poly{Int}) || push!(fails, "p - p $trial")
    x = vals()
    note("evaluate(p*q) vs evaluate(p)*evaluate(q)",
         relerr(S.evaluate(p * q, x), S.evaluate(p, x) * S.evaluate(q, x)))
    # substitution by power products with rational exponents
    sub = Dict(v => S.mono(rand(rng, cons) => randexp(false), keq => randexp(false))
               for v in cons if rand(rng, Bool))
    y = copy(x)
    for (v, mm) in sub
        y[v] = prod((x[w]^BigFloat(e) for (w, e) in mm); init = big(1))
    end
    note("subst", relerr(S.evaluate(S.subst(p, v -> get(sub, v, nothing)), x),
                         S.evaluate(p, y)))
    # exact division by a Haldane-like binomial, quotient recovered and remainder caught
    d = S.term(Int, S.mono(conc[1] => 1, conc[2] => 1)) -
        S.term(Int, S.mono(conc[3] => 1, conc[4] => 1, keq => -1 // rand(rng, 1:2)))
    lq = randpoly(Int; laurent = true)
    iszero(lq) && continue
    S.divexact(lq * d, d) == lq || push!(fails, "divexact quotient $trial")
    bad = lq * d + S.term(Int, S.mono(cons[1] => 1))
    try
        S.divexact(bad, d)
        push!(fails, "divexact accepted a remainder $trial")
    catch e
        e isa ArgumentError || push!(fails, "divexact error type $trial: $e")
    end
    # lowest terms: the shift zeroes every concentration's lowest exponent
    num, den = randpoly(Int; laurent = true), randpoly(Int; laurent = true)
    (iszero(num) || iszero(den)) && continue
    sh = S.conc_shift((num, den), vars)
    n2, d2 = num * sh, den * sh
    for v in conc
        lo = minimum(S.mexp(m, v) for p in (n2, d2) for m in keys(p.terms))
        lo == 0 || push!(fails, "conc_shift leaves exponent $lo $trial")
    end
    note("conc_shift ratio", relerr(S.evaluate(n2, x) / S.evaluate(d2, x),
                                    S.evaluate(num, x) / S.evaluate(den, x)))
    # coefficients by concentration monomial recombine to the polynomial
    parts = S.split_conc(p, vars)
    back = reduce(+, (c * cm for (cm, c) in parts); init = S.Poly{Int}())
    back == p || push!(fails, "split_conc $trial")
end
# Float64 evaluation agrees with BigFloat
p = randpoly(Int; n = 6)
x = vals()
note("Float64 vs BigFloat", relerr(big(S.evaluate(p, Float64.(x))), S.evaluate(p, x)))
# checked Int arithmetic overflows; BigInt takes over
big_p = S.term(Int, S.mono(cons[1] => 1), typemax(Int) ÷ 2) +
        S.term(Int, S.mono(cons[2] => 1), 3)
overflowed = try
    big_p * big_p
    false
catch e
    e isa OverflowError
end
widened = convert(S.Poly{BigInt}, big_p) * convert(S.Poly{BigInt}, big_p)
@printf("overflow: Int product throws OverflowError %s; BigInt product has %d terms\n",
        overflowed, S.nterms(widened))
println("failures: ", length(fails))
foreach(f -> println("    ", f), first(fails, 10))
for (k, v) in sort!(collect(worst))
    @printf("worst relative error, %-45s %.1e\n", k, Float64(v))
end
