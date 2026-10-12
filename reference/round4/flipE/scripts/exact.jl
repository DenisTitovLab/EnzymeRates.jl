# ABOUTME: Task 2: flip-at-free-E rate laws of iso mechanisms (a second free form Eprime or D)
# ABOUTME: against the two-conformation mass-action oracle in which only free E flips.
include(joinpath(@__DIR__, "lib.jl"))

ISO_UU_SS(tS, tC, tP, tI) = [("E + S <--> E(S)", tS), ("E(S) <--> Eprime(P)", tC),
                                   ("Eprime + P <--> Eprime(P)", tP), ("Eprime <--> E", tI)]
ISO_UU_RE(tS, tC, tP, tI) = [("E + S ⇌ E(S)", tS), ("E(S) <--> Eprime(P)", tC),
                                   ("Eprime + P ⇌ Eprime(P)", tP), ("Eprime <--> E", tI)]
ISO_BIBI(tA, tB, tC, tP, tQ, tI) = [
    ("E + A ⇌ E(A)", tA), ("E(A) + B <--> E(A, B)", tB),
    ("E(A, B) <--> Eprime(P, Q)", tC), ("Eprime(Q) + P <--> Eprime(P, Q)", tP),
    ("Eprime + Q <--> Eprime(Q)", tQ), ("Eprime <--> E", tI)]
const EQ, NE, OA = "EqualAI", "NonequalAI", "OnlyA"

cases = Tuple{String, Vector{String}, Vector{String}, Vector{Tuple{String, String}}, Int}[]
for (lbl, f) in (("IsoUU-SS", ISO_UU_SS), ("IsoUU-RE", ISO_UU_RE))
    for (name, tags) in (("iso NonequalAI", (EQ, EQ, EQ, NE)),
                         ("catalysis NonequalAI", (EQ, NE, EQ, EQ)),
                         ("catalysis OnlyA", (EQ, OA, EQ, EQ)),
                         ("S binding NonequalAI", (NE, EQ, EQ, EQ)),
                         ("iso OnlyA", (EQ, EQ, EQ, OA)),
                         ("catalysis OnlyA + iso NonequalAI", (EQ, OA, EQ, NE)),
                         ("all NonequalAI", (NE, NE, NE, NE)),
                         ("S binding + catalysis OnlyA", (OA, OA, EQ, EQ)))
        push!(cases, ("$lbl: $name", ["S"], ["P"], f(tags...), 1))
    end
end
for (name, tags) in (("iso NonequalAI", (EQ, EQ, EQ, EQ, EQ, NE)),
                     ("catalysis NonequalAI", (EQ, EQ, NE, EQ, EQ, EQ)),
                     ("catalysis OnlyA", (EQ, EQ, OA, EQ, EQ, EQ)),
                     ("B binding NonequalAI", (EQ, NE, EQ, EQ, EQ, EQ)),
                     ("iso OnlyA", (EQ, EQ, EQ, EQ, EQ, OA)),
                     ("A binding (RE) NonequalAI", (NE, EQ, EQ, EQ, EQ, EQ)),
                     ("P release + iso NonequalAI", (EQ, EQ, EQ, NE, EQ, NE)),
                     ("A binding + catalysis OnlyA", (OA, EQ, OA, EQ, EQ, EQ)),
                     ("all NonequalAI", (NE, NE, NE, NE, NE, NE)))
    push!(cases, ("IsoBiBi: $name", ["A", "B"], ["P", "Q"], ISO_BIBI(tags...), 1))
end
# A second free form D in rapid equilibrium with E sorts first and roots E's segment, so
# free E's weight carries its alpha.
push!(cases, ("IsoUU-D: D ⇌ E NonequalAI", ["S"], ["P"],
              [ISO_UU_SS(EQ, EQ, EQ, EQ); ("D ⇌ E", NE)], 1))
push!(cases, ("IsoUU-D: D ⇌ E NonequalAI, catalysis NonequalAI", ["S"], ["P"],
              [ISO_UU_SS(EQ, NE, EQ, EQ); ("D ⇌ E", NE)], 1))
# Two catalytic subunits (concerted oracle over protomer pairs).
push!(cases, ("IsoUU-SS n=2: iso NonequalAI", ["S"], ["P"], ISO_UU_SS(EQ, EQ, EQ, NE), 2))
push!(cases, ("IsoUU-SS n=2: catalysis OnlyA + iso NonequalAI", ["S"], ["P"],
              ISO_UU_SS(EQ, OA, EQ, NE), 2))
push!(cases, ("IsoBiBi n=2: P release + iso NonequalAI", ["A", "B"], ["P", "Q"],
              ISO_BIBI(EQ, EQ, EQ, NE, EQ, NE), 2))

println("package: ", pathof(EnzymeRates))
@printf("%-52s %9s %9s %9s %9s %9s %s\n", "case", "worst", "eq-self", "L-sens",
        "alt-flip", "mass_I", "notes")
summary = Tuple{String, BigFloat}[]
for (label, subs, prods, steps, n) in cases
    m = build(dsl(subs, prods, steps; n))
    if m isa String
        println(rpad(label, 52), " NOT CONSTRUCTED: ", m)
        continue
    end
    am = ER.AllostericMechanism(m)
    fp = ER.fitted_params(m)
    x1 = Symbol(subs[1])
    rng = MersenneTwister(hash(label) % UInt32)
    worst, eqself, lsens, flipdiff, massI = big(0), big(0), big(Inf), big(0), big(0)
    notes = String[]
    pruned_pkg = inactive_forms(am)
    for draw in 1:(n == 1 ? 5 : 2)
        fitted = draw_fitted(rng, fp)
        Keq = exp(BigFloat(log(0.5) + log(40) * rand(rng)))
        vals = all_values(am, fitted, Keq)
        sA, sI, miss = oracle_steps(am, vals)
        isempty(miss) || push!(notes, "missing I constants: $miss")
        for point in 1:3
            concs = Dict{Symbol, BigFloat}(Symbol(x) => draw_conc(rng) for x in [subs; prods])
            v = pkg_rate(m, vals, concs)
            vo, left, mI = oracle_v(sA, sI, concs, vals[:L], x1; n)
            worst = max(worst, relerr(v, vo))
            massI = max(massI, mI)
            # The oracle's inactive forms must be the package's pruned ones.
            keptI = sort!(collect(setdiff(Set(f for e in sI for f in (e.a, e.b)), left)))
            keptI == pruned_pkg || push!(notes, "inactive forms differ: oracle $keptI, package $pruned_pkg")
            # Does L matter here, and would a flip at Eprime be told apart?
            v0, _, _ = oracle_v(sA, sI, concs, big(0), x1; n)
            lsens = min(lsens, abs(vo - v0) / abs(vo))
            # The law against an oracle flipping at the other free form (D when present,
            # else Eprime): large values show the oracle tells the flip form apart.
            alt = any(e -> :D in (e.a, e.b), sA) ? :D : :Eprime
            if n == 1
                vp, _, _ = oracle_v(sA, sI, concs, vals[:L], x1; n, flip = alt)
                flipdiff = max(flipdiff, abs(v - vp) / abs(vo))
            end
            # Self-check of the oracle: no net rate at equilibrium.
            ce = eq_concs(rng, am, Keq, concs)
            ve, _, _ = oracle_v(sA, sI, ce, vals[:L], x1; n)
            eqself = max(eqself, abs(ve) / abs(vo))
        end
    end
    push!(summary, (label, worst))
    @printf("%-52s %9.1e %9.1e %9.1e %9.1e %9.2f %s\n", label, worst, eqself, lsens,
            flipdiff, massI, join(unique(notes), "; "))
end
println("\nworst over all cases: ", @sprintf("%.2e", maximum(last.(summary))),
        " over ", length(summary), " constructed cases")
