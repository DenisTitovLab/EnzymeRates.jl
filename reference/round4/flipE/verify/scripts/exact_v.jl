# ABOUTME: Checks the checkout's allosteric laws for mechanisms with a second free form against
# ABOUTME: the independent oracle (vlib.jl) flipping at E and at each other free form.
include(joinpath(@__DIR__, "vlib.jl"))

const DF = Dict("ES" => "E(S)", "EP" => "E(P)", "EprimeP" => "Eprime(P)", "EA" => "E(A)",
                "EAB" => "E(A, B)", "EprimePQ" => "Eprime(P, Q)", "EprimeQ" => "Eprime(Q)",
                "FP" => "F(P)")
V(f, c, t, r, re, tag) = VStep(f, c, t, r, re, tag)
const EQ, NE, OA = :EqualAI, :NonequalAI, :OnlyA

uu_ss(tS, tC, tP, tI) = [V("E", ["S"], "ES", String[], false, tS),
                         V("ES", String[], "EprimeP", String[], false, tC),
                         V("Eprime", ["P"], "EprimeP", String[], false, tP),
                         V("Eprime", String[], "E", String[], false, tI)]
uu_re(tS, tC, tP, tI) = [V("E", ["S"], "ES", String[], true, tS),
                         V("ES", String[], "EprimeP", String[], false, tC),
                         V("Eprime", ["P"], "EprimeP", String[], true, tP),
                         V("Eprime", String[], "E", String[], false, tI)]
bibi(tA, tB, tC, tP, tQ, tI) = [V("E", ["A"], "EA", String[], true, tA),
                                V("EA", ["B"], "EAB", String[], false, tB),
                                V("EAB", String[], "EprimePQ", String[], false, tC),
                                V("EprimeQ", ["P"], "EprimePQ", String[], false, tP),
                                V("Eprime", ["Q"], "EprimeQ", String[], false, tQ),
                                V("Eprime", String[], "E", String[], false, tI)]
# A plain uni-uni on E with a dead-end free form X joined to E (RE or SS).
uu_dead(X, re, tX; tC = EQ) = [V("E", ["S"], "ES", String[], false, EQ),
                               V("ES", String[], "EP", String[], false, tC),
                               V("E", ["P"], "EP", String[], false, EQ),
                               V(X, String[], "E", String[], re, tX)]
# Rapid-equilibrium random uni-uni (one segment) with a free form Ca ⇌ E sorting first.
uu_one_segment = [V("E", ["S"], "ES", String[], true, EQ),
                  V("ES", String[], "EP", String[], false, NE),
                  V("E", ["P"], "EP", String[], true, EQ),
                  V("Ca", String[], "E", String[], true, NE)]
# Three free forms on one cycle: E → ES → F(P) → F → G → E.
three_free(tC, tFG, tGE) = [V("E", ["S"], "ES", String[], true, EQ),
                            V("ES", String[], "FP", String[], false, tC),
                            V("F", ["P"], "FP", String[], true, EQ),
                            V("F", String[], "G", String[], false, tFG),
                            V("G", String[], "E", String[], false, tGE)]

cases = [
    ("IsoUU-SS iso NE", ["S"], ["P"], uu_ss(EQ, EQ, EQ, NE), [1, 2, 3]),
    ("IsoUU-SS cat OA + iso NE", ["S"], ["P"], uu_ss(EQ, OA, EQ, NE), [1, 2]),
    ("IsoUU-RE all NE", ["S"], ["P"], uu_re(NE, NE, NE, NE), [1]),
    ("IsoUU-SS iso OA", ["S"], ["P"], uu_ss(EQ, EQ, EQ, OA), [1, 2]),
    ("IsoBiBi P-release NE + iso NE", ["A", "B"], ["P", "Q"], bibi(EQ, EQ, EQ, NE, EQ, NE), [1, 2]),
    ("IsoBiBi cat OA", ["A", "B"], ["P", "Q"], bibi(EQ, EQ, OA, EQ, EQ, EQ), [1]),
    ("Dead-end D ⇌ E NE", ["S"], ["P"], uu_dead("D", true, NE), [1, 2]),
    ("Dead-end D ⇌ E NE, cat NE", ["S"], ["P"], uu_dead("D", true, NE; tC = NE), [1]),
    ("Dead-end D <--> E (SS) NE", ["S"], ["P"], uu_dead("D", false, NE), [1]),
    ("One RE segment, Ca ⇌ E NE", ["S"], ["P"], uu_one_segment, [1, 2]),
    ("Three free forms, F<-->G NE", ["S"], ["P"], three_free(EQ, NE, EQ), [1, 2]),
    ("Three free forms, cat NE, G<-->E OA", ["S"], ["P"], three_free(NE, EQ, OA), [1]),
]

println("package: ", pathof(EnzymeRates))
@printf("%-38s %2s %10s %10s %10s %-22s %s\n", "case", "n", "vs flip@E", "eq-self", "L-sens",
        "vs flip@other", "notes")
for (label, subs, prods, steps, ns) in cases, n in ns
    m = build(dsltext(subs, prods, steps, DF; n))
    if m isa String
        @printf("%-38s %2d NOT CONSTRUCTED: %s\n", label, n, first(replace(m, "\n" => " "), 140))
        continue
    end
    fp = String.(ER.fitted_params(m))
    rng = MersenneTwister(hash(label) % UInt32 + n)
    others = filter(!=("E"), [f for f in ("Eprime", "D", "Ca", "F", "G")
                              if any(s -> f in (s.from, s.to), steps)])
    worstE, eqself, lsens = big(0.0), big(0.0), big(Inf)
    worstO = Dict{String, BigFloat}(o => big(0.0) for o in others)
    notes = String[]
    for draw in 1:(n == 1 ? 3 : 2)
        fitted = Dict(p => (p == "L" ? lu(rng, 0.2, 5) : lu(rng, 0.1, 10)) for p in fp)
        Keq = lu(rng, 0.5, 20)
        vals = try
            solve_values(steps, fitted, Keq, subs, prods)
        catch err
            push!(notes, "solve: " * first(sprint(showerror, err), 120)); break
        end
        L = fitted["L"]
        pts = Dict{String, BigFloat}[]
        for k in 1:2
            push!(pts, Dict(x => lu(rng, 0.1, 10) for x in [subs; prods]))
        end
        push!(pts, merge(Dict(x => lu(rng, 0.1, 10) for x in subs),
                         Dict(x => big(0) for x in prods)))
        push!(pts, merge(Dict(x => big(10)^6 for x in subs),
                         Dict(x => lu(rng, 0.1, 10) for x in prods)))
        # Errors are relative to the larger of the two rates and the generic point's rate,
        # so a zero rate against roundoff does not count as a mismatch.
        typ = abs(oracle_rate(steps, vals, pts[1], L, subs[1]; n)[1])
        for c in pts
            v = pkg_rate(m, fitted, Keq, c)
            vo, _ = oracle_rate(steps, vals, c, L, subs[1]; n)
            scale = max(abs(vo), abs(v), typ)
            err = abs(v - vo) / scale
            worstE = max(worstE, err)
            v0, _ = oracle_rate(steps, vals, c, big(0), subs[1]; n)
            vo != 0 && (lsens = min(lsens, abs(vo - v0) / abs(vo)))
            for o in others
                vx, _ = oracle_rate(steps, vals, c, L, subs[1]; n, flip = o)
                worstO[o] = max(worstO[o], abs(v - vx) / max(abs(vo), abs(vx), typ))
            end
            all(iszero(c[x]) for x in prods) && vo == 0 && push!(notes, "v=0 at products 0")
        end
        # Equilibrium self-check of the oracle: last product set so ΠP = Keq·ΠS.
        ce = Dict(x => lu(rng, 0.1, 10) for x in [subs; prods])
        ce[prods[end]] = Keq * prod(ce[x] for x in subs) /
                         prod((ce[x] for x in prods[1:end-1]); init = big(1))
        ve, _ = oracle_rate(steps, vals, ce, L, subs[1]; n)
        vt, _ = oracle_rate(steps, vals, pts[1], L, subs[1]; n)
        eqself = max(eqself, abs(ve) / abs(vt))
    end
    o = join(["$k:" * @sprintf("%.1e", Float64(v)) for (k, v) in worstO], " ")
    @printf("%-38s %2d %10.1e %10.1e %10.1e %-22s %s\n", label, n, Float64(worstE),
            Float64(eqself), Float64(lsens), o, join(unique(notes), "; "))
end
