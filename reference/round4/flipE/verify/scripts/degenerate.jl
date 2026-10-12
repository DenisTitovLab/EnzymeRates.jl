# ABOUTME: Checks whether mechanisms newly admitted by flip-at-E (an :OnlyA step inside free
# ABOUTME: E's cycle) give a law that vanishes at zero products, and what _kcat_forward does.
using EnzymeRates, Random, Printf
const ER = EnzymeRates

mk(src) = try
    Core.eval(Main, Meta.parse(src))
catch err
    err isa LoadError && (err = err.error)
    sprint(showerror, err)
end

dsl(subs, prods, steps; n = 1) =
    "@allosteric_mechanism begin\n    substrates: $(join(subs, ", "))\n" *
    "    products: $(join(prods, ", "))\n    catalytic_multiplicity: $n\n" *
    "    catalytic_steps: begin\n" * join(["        $s :: $t" for (s, t) in steps], "\n") *
    "\n    end\nend"

cases = [
    ("UU-SS iso OnlyA", ["S"], ["P"],
     [("E + S <--> E(S)", "EqualAI"), ("E(S) <--> Eprime(P)", "EqualAI"),
      ("Eprime + P <--> Eprime(P)", "EqualAI"), ("Eprime <--> E", "OnlyA")]),
    ("UU-RE iso OnlyA", ["S"], ["P"],
     [("E + S ⇌ E(S)", "EqualAI"), ("E(S) <--> Eprime(P)", "EqualAI"),
      ("Eprime + P ⇌ Eprime(P)", "EqualAI"), ("Eprime <--> E", "OnlyA")]),
    ("UU-SS iso OnlyA, catalysis NonequalAI", ["S"], ["P"],
     [("E + S <--> E(S)", "EqualAI"), ("E(S) <--> Eprime(P)", "NonequalAI"),
      ("Eprime + P <--> Eprime(P)", "EqualAI"), ("Eprime <--> E", "OnlyA")]),
    ("BiBi iso OnlyA", ["A", "B"], ["P", "Q"],
     [("E + A ⇌ E(A)", "EqualAI"), ("E(A) + B <--> E(A, B)", "EqualAI"),
      ("E(A, B) <--> Eprime(P, Q)", "EqualAI"), ("Eprime(Q) + P <--> Eprime(P, Q)", "EqualAI"),
      ("Eprime + Q <--> Eprime(Q)", "EqualAI"), ("Eprime <--> E", "OnlyA")]),
    # Controls: the same shapes with the isomerization kept in both conformations.
    ("UU-SS iso NonequalAI (control)", ["S"], ["P"],
     [("E + S <--> E(S)", "EqualAI"), ("E(S) <--> Eprime(P)", "EqualAI"),
      ("Eprime + P <--> Eprime(P)", "EqualAI"), ("Eprime <--> E", "NonequalAI")]),
    ("UU + D ⇌ E NonequalAI (alpha case)", ["S"], ["P"],
     [("E + S <--> E(S)", "EqualAI"), ("E(S) <--> E(P)", "EqualAI"),
      ("E + P <--> E(P)", "EqualAI"), ("D ⇌ E", "NonequalAI")]),
]

println("package: ", pathof(EnzymeRates))
for (label, subs, prods, steps) in cases, n in (1, 2)
    m = mk(dsl(subs, prods, steps; n))
    println("\n=== $label  n=$n")
    if m isa String
        println("  NOT CONSTRUCTED: ", first(m, 300))
        continue
    end
    em = m
    fp = ER.fitted_params(em)
    n == 1 && println("  fitted_params: ", fp)
    mets = Tuple(ER.metabolites(em))
    rng = MersenneTwister(7)
    for draw in 1:3
        vals = Tuple(p === :L ? exp(log(0.2) + log(25) * rand(rng)) :
                     exp(log(0.1) + log(100) * rand(rng)) for p in fp)
        prm = NamedTuple{(fp..., :Keq, :E_total)}((vals..., 5.0, 1.0))
        c(f) = NamedTuple{mets}(Tuple(f(x) for x in mets))
        isS(x) = String(x) in subs
        v0 = ER.rate_equation(em, c(x -> isS(x) ? 2.0 : 0.0), prm)
        vsmall = ER.rate_equation(em, c(x -> isS(x) ? 2.0 : 1e-6), prm)
        v1 = ER.rate_equation(em, c(x -> isS(x) ? 2.0 : 1.0), prm)
        vsat = ER.rate_equation(em, c(x -> isS(x) ? 1e9 : 0.0), prm)
        vsatP = ER.rate_equation(em, c(x -> isS(x) ? 1e9 : 1e-6), prm)
        kc = try
            @sprintf("%.6g", ER._kcat_forward(em, prm))
        catch err
            "ERROR: " * first(sprint(showerror, err), 160)
        end
        @printf("  draw %d: v(S=2, products 0) = %.6g | products 1e-6: %.6g | products 1: %.6g\n",
                draw, v0, vsmall, v1)
        @printf("          v(S=1e9, products 0) = %.6g | v(S=1e9, products 1e-6) = %.6g | _kcat_forward = %s\n",
                vsat, vsatP, kc)
    end
end
