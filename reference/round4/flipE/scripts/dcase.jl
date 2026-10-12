# ABOUTME: A second free form D ⇌ E (RE, :NonequalAI) roots E's segment: which form the law
# ABOUTME: flips at, against oracles flipping at E and at D (run under both checkouts).
include(joinpath(@__DIR__, "lib.jl"))

cases = [
    ("UU + D ⇌ E NonequalAI", [("E + S <--> E(S)", "EqualAI"), ("E(S) <--> E(P)", "EqualAI"),
                               ("E + P <--> E(P)", "EqualAI"), ("D ⇌ E", "NonequalAI")]),
    ("UU + D ⇌ E NonequalAI, catalysis NonequalAI",
     [("E + S <--> E(S)", "EqualAI"), ("E(S) <--> E(P)", "NonequalAI"),
      ("E + P <--> E(P)", "EqualAI"), ("D ⇌ E", "NonequalAI")]),
    ("UU + Ez ⇌ E NonequalAI (Ez sorts after E)",
     [("E + S <--> E(S)", "EqualAI"), ("E(S) <--> E(P)", "EqualAI"),
      ("E + P <--> E(P)", "EqualAI"), ("Ez ⇌ E", "NonequalAI")]),
]
println("package: ", pathof(EnzymeRates))
for (label, steps) in cases
    m = build(dsl(["S"], ["P"], steps))
    m isa String && (println(label, ": NOT CONSTRUCTED: ", m); continue)
    am = ER.AllostericMechanism(m)
    other = first(f for f in (:D, :Ez) if any(s -> ER.name(ER.from_species(s)) == f ||
                                                 ER.name(ER.to_species(s)) == f,
                                              Iterators.flatten(ER.steps(am))))
    rng = MersenneTwister(11)
    eE, eO = big(0), big(0)
    for draw in 1:4
        vals = all_values(am, draw_fitted(rng, ER.fitted_params(m)),
                          exp(BigFloat(log(0.5) + log(40) * rand(rng))))
        sA, sI, _ = oracle_steps(am, vals)
        for _ in 1:3
            concs = Dict(:S => draw_conc(rng), :P => draw_conc(rng))
            v = pkg_rate(m, vals, concs)
            vE, _, _ = oracle_v(sA, sI, concs, vals[:L], :S)
            vO, _, _ = oracle_v(sA, sI, concs, vals[:L], :S; flip = other)
            eE, eO = max(eE, relerr(v, vE)), max(eO, relerr(v, vO))
        end
    end
    @printf("%-46s law vs flip-at-E oracle %.1e | vs flip-at-%s oracle %.1e\n", label, eE,
            other, eO)
end
