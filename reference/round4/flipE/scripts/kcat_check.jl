# ABOUTME: kcat of the newly constructible two-free-form mechanisms: the package's
# ABOUTME: _kcat_forward against its own rate at saturating substrates with products at 0.
include(joinpath(@__DIR__, "lib.jl"))
include(joinpath(@__DIR__, "fixtures_def.jl"))

EQ, NE, OA = "EqualAI", "NonequalAI", "OnlyA"
uu(tS, tC, tP, tI) = [("E + S <--> E(S)", tS), ("E(S) <--> Eprime(P)", tC),
                      ("Eprime + P <--> Eprime(P)", tP), ("Eprime <--> E", tI)]
bibi(tC, tI) = [("E + A ⇌ E(A)", EQ), ("E(A) + B <--> E(A, B)", EQ),
                ("E(A, B) <--> Eprime(P, Q)", tC), ("Eprime(Q) + P <--> Eprime(P, Q)", EQ),
                ("Eprime + Q <--> Eprime(Q)", EQ), ("Eprime <--> E", tI)]
cases = vcat([("IsoUU iso NonequalAI", ["S"], ["P"], String[], uu(EQ, EQ, EQ, NE)),
              ("IsoUU catalysis OnlyA", ["S"], ["P"], String[], uu(EQ, OA, EQ, EQ)),
              ("IsoUU catalysis OnlyA + iso NonequalAI", ["S"], ["P"], String[],
               uu(EQ, OA, EQ, NE)),
              ("IsoUU all NonequalAI", ["S"], ["P"], String[], uu(NE, NE, NE, NE)),
              ("IsoBiBi iso NonequalAI", ["A", "B"], ["P", "Q"], String[], bibi(EQ, NE)),
              ("IsoBiBi catalysis OnlyA", ["A", "B"], ["P", "Q"], String[], bibi(OA, EQ))],
             fixtures)
println("package: ", pathof(EnzymeRates))
for (label, subs, prods, inh, steps) in cases, n in (1, 2)
    m = build(dsl(subs, prods, steps; n, inhibitors = inh))
    m isa String && (println(rpad(label, 44), " n=$n NOT CONSTRUCTED: ", m); continue)
    rng = MersenneTwister(5)
    worst, kc1 = 0.0, NaN
    for _ in 1:4
        fp = ER.fitted_params(m)
        prm = NamedTuple{(fp..., :Keq, :E_total)}(
            (Tuple(p === :L ? exp(log(0.2) + log(25) * rand(rng)) :
                   exp(log(0.1) + log(100) * rand(rng)) for p in fp)..., 3.0, 1.0))
        kc = ER._kcat_forward(m, prm)
        mets = Tuple(ER.metabolites(m))
        sat = NamedTuple{mets}(Tuple(Symbol(x) in Symbol.(subs) ? 1e9 : 0.0 for x in String.(mets)))
        v = ER.rate_equation(m, sat, prm)
        worst = max(worst, abs(kc - v) / abs(v))
        kc1 = kc
    end
    @printf("%-44s n=%d  max |kcat - v(substrates 1e9, products 0)|/v = %.1e  (kcat %.4g)\n",
            label, n, worst, kc1)
end
