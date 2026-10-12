# ABOUTME: For free forms joined to E in rapid equilibrium (D ⇌ E, Ca ⇌ E), prints the law at L and
# ABOUTME: at L·(K_I/K_A)^n, to test whether flip-at-D and flip-at-E laws differ only by L's scale.
include(joinpath(@__DIR__, "vlib.jl"))

const DF = Dict("ES" => "E(S)", "EP" => "E(P)")
V(f, c, t, r, re, tag) = VStep(f, c, t, r, re, tag)
const EQ, NE = :EqualAI, :NonequalAI
dead(X) = [V("E", ["S"], "ES", String[], false, EQ), V("ES", String[], "EP", String[], false, NE),
           V("E", ["P"], "EP", String[], false, EQ), V(X, String[], "E", String[], true, NE)]
seg(X) = [V("E", ["S"], "ES", String[], true, EQ), V("ES", String[], "EP", String[], false, NE),
          V("E", ["P"], "EP", String[], true, EQ), V(X, String[], "E", String[], true, NE)]

println("package: ", pathof(EnzymeRates))
for (label, steps, X) in (("dead-end D ⇌ E", dead("D"), "D"), ("one segment Ca ⇌ E", seg("Ca"), "Ca")),
    n in (1, 2)
    m = build(dsltext(["S"], ["P"], steps, DF; n))
    m isa String && (println(label, " NOT CONSTRUCTED"); continue)
    fp = ER.fitted_params(m)
    rng = MersenneTwister(42 + n)
    for draw in 1:2
        vals = Dict(String(p) => (p === :L ? lu(rng, 0.2, 5) : lu(rng, 0.1, 10)) for p in fp)
        Keq = lu(rng, 0.5, 20)
        r = (vals["K_I_$(X)_to_E"] / vals["K_A_$(X)_to_E"])^n
        for (S, P) in ((big(0.7), big(1.3)), (big(4.0), big(0.0)))
            c = Dict("S" => S, "P" => P)
            v1 = pkg_rate(m, vals, Keq, c)
            vals2 = copy(vals); vals2["L"] = vals["L"] * r
            v2 = pkg_rate(m, vals2, Keq, c)
            @printf("%-20s n=%d draw %d S=%.1f P=%.1f  law(L) %.15e  law(L·(K_I/K_A)^n) %.15e\n",
                    label, n, draw, S, P, v1, v2)
        end
    end
end
