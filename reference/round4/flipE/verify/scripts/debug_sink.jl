# ABOUTME: Point-by-point law vs oracle for the two cases where exact_v.jl reports a mismatch
# ABOUTME: (an :OnlyA step breaking the inactive cycle, leaving a product-dependent sink).
include(joinpath(@__DIR__, "vlib.jl"))

const DF = Dict("ES" => "E(S)", "EprimeP" => "Eprime(P)", "FP" => "F(P)")
V(f, c, t, r, re, tag) = VStep(f, c, t, r, re, tag)
const EQ, NE, OA = :EqualAI, :NonequalAI, :OnlyA
uu_ss(tS, tC, tP, tI) = [V("E", ["S"], "ES", String[], false, tS),
                         V("ES", String[], "EprimeP", String[], false, tC),
                         V("Eprime", ["P"], "EprimeP", String[], false, tP),
                         V("Eprime", String[], "E", String[], false, tI)]
three_free(tC, tFG, tGE) = [V("E", ["S"], "ES", String[], true, EQ),
                            V("ES", String[], "FP", String[], false, tC),
                            V("F", ["P"], "FP", String[], true, EQ),
                            V("F", String[], "G", String[], false, tFG),
                            V("G", String[], "E", String[], false, tGE)]

for (label, steps, n) in (("IsoUU-SS iso OA", uu_ss(EQ, EQ, EQ, OA), 2),
                          ("IsoUU-SS iso OA", uu_ss(EQ, EQ, EQ, OA), 1),
                          ("Three free forms, cat NE, G<-->E OA", three_free(NE, EQ, OA), 1),
                          ("Three free forms, cat NE, G<-->E OA", three_free(NE, EQ, OA), 2))
    m = build(dsltext(["S"], ["P"], steps, DF; n))
    println("\n=== $label n=$n")
    m isa String && (println("NOT CONSTRUCTED: ", m); continue)
    fp = String.(ER.fitted_params(m))
    rng = MersenneTwister(3)
    fitted = Dict(p => (p == "L" ? lu(rng, 0.2, 5) : lu(rng, 0.1, 10)) for p in fp)
    Keq = lu(rng, 0.5, 20)
    vals = solve_values(steps, fitted, Keq, ["S"], ["P"])
    for (S, P) in ((2.0, 3.0), (2.0, 0.5), (2.0, 1e-3), (2.0, 0.0), (1e6, 3.0), (1e6, 0.5))
        c = Dict("S" => big(S), "P" => big(P))
        v = pkg_rate(m, fitted, Keq, c)
        vo, mI = oracle_rate(steps, vals, c, fitted["L"], "S"; n)
        @printf("  S=%-8g P=%-6g law %+.10e  oracle %+.10e  rel %.1e  oracle I-mass %.6f\n",
                S, P, Float64(v), Float64(vo), Float64(abs(v - vo) / max(abs(vo), abs(v), big(10)^-300)),
                Float64(mI))
    end
end
