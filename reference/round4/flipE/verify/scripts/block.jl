# ABOUTME: R4 (dead-end F ⇌ E :: OnlyA) and R2 (OnlyA binding + OnlyA isomerization) with the
# ABOUTME: :OnlyA steps kept in the inactive state but slowed by ε, against the deletion law.
include(joinpath(@__DIR__, "vlib.jl"))

const DF = Dict("ES" => "E(S)", "EP" => "E(P)")
V(f, c, t, r, re, tag) = VStep(f, c, t, r, re, tag)
const EQ, NE, OA = :EqualAI, :NonequalAI, :OnlyA

r4(tF) = [V("E", ["S"], "ES", String[], true, EQ), V("ES", String[], "EP", String[], false, NE),
          V("E", ["P"], "EP", String[], true, EQ), V("F", String[], "E", String[], true, tF)]
r2 = [V("E", ["S"], "ES", String[], true, OA), V("ES", String[], "EP", String[], false, EQ),
      V("F", ["P"], "EP", String[], true, EQ), V("F", String[], "E", String[], false, OA)]

println("package: ", pathof(EnzymeRates))
println("\n## R4: F ⇌ E :: OnlyA. Kinetic block keeps F ⇌ E in the inactive state at rates ε·(ratio, 1).")
for n in (1, 2)
    m4 = build(dsltext(["S"], ["P"], r4(OA), DF; n))
    m5 = build(dsltext(["S"], ["P"], r4(NE), DF; n))
    rng = MersenneTwister(100 + n)
    for draw in 1:2
        fitted = Dict(String(p) => (p === :L ? lu(rng, 0.2, 5) : lu(rng, 0.1, 10))
                      for p in ER.fitted_params(m5))
        Keq = lu(rng, 0.5, 20)
        KA = haskey(fitted, "K_A_E_to_F") ? "K_A_E_to_F" : "K_A_F_to_E"
        f4 = Dict(k => v for (k, v) in fitted if String(k) in String.(ER.fitted_params(m4)))
        f5 = copy(fitted); f5[replace(KA, "K_A_" => "K_I_")] = fitted[KA]   # EqualAI twin
        vals = solve_values(r4(OA), f4, Keq, ["S"], ["P"])
        ratio = stepconsts(r4(OA)[4], :A, y -> haskey(vals, y)).ratio(y -> vals[y])
        for (S, P) in ((big(0.8), big(1.5)), (big(3.0), big(0.0)))
            c = Dict("S" => S, "P" => P)
            vdel = pkg_rate(m4, f4, Keq, c)
            vtwin = pkg_rate(m5, f5, Keq, c)
            line = @sprintf("n=%d draw %d S=%.1f P=%.1f  deletion law %.10f  EqualAI-twin law %.10f |",
                            n, draw, S, P, vdel, vtwin)
            for ε in (big(1), big(10)^-6, big(10)^-14)
                vb, _ = oracle_rate(r4(OA), vals, c, fitted["L"], "S"; n,
                                    blocked = Dict(4 => (ε * ratio, ε)))
                line *= @sprintf(" ε=%.0e: vs twin %.1e vs deletion %.1e |", ε,
                                 abs(vb - vtwin) / abs(vtwin), abs(vb - vdel) / abs(vdel))
            end
            println(line)
        end
    end
end

println("\n## R2: E + S ⇌ E(S) :: OnlyA, F <--> E :: OnlyA. Kinetic block: the binding keeps")
println("## rapid equilibrium with K/ε; the isomerization keeps both rates × ε^a.")
for n in (1,)
    m2 = build(dsltext(["S"], ["P"], r2, DF; n))
    m2 isa String && (println("R2 NOT CONSTRUCTED: ", m2); continue)
    rng = MersenneTwister(7)
    for draw in 1:2
        fitted = Dict(String(p) => (p === :L ? lu(rng, 0.2, 5) : lu(rng, 0.1, 10))
                      for p in ER.fitted_params(m2))
        Keq = lu(rng, 0.5, 20)
        vals = solve_values(r2, fitted, Keq, ["S"], ["P"])
        K = vals["K_A_ES_to_E_S"]
        kFE, kEF = vals["k_A_F_to_E"], vals["k_A_E_to_F"]
        for (S, P) in ((big(2.0), big(1.0)), (big(2.0), big(10)^-6), (big(2.0), big(0.0)))
            c = Dict("S" => S, "P" => P)
            vdel = pkg_rate(m2, fitted, Keq, c)
            line = @sprintf("draw %d S=%.1f P=%.0e deletion law %.6e |", draw, S, P, vdel)
            for a in (0.5, 1, 2), ε in (big(10)^-6, big(10)^-12)
                Λb = big(10)^30
                vb, _ = oracle_rate(r2, vals, c, fitted["L"], "S"; n,
                                    blocked = Dict(1 => (Λb * ε / K, Λb), 4 => (ε^a * kFE, ε^a * kEF)))
                line *= @sprintf(" a=%.1f ε=%.0e: %.2e", a, ε, abs(vb - vdel) / abs(vb))
            end
            println(line)
        end
    end
end
