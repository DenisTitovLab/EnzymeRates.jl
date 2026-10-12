# ABOUTME: Checks the guard probe's RE-turnover, face-class, D=0-face and N≡0 claims with an
# ABOUTME: independent oracle on free-energy constants (stdlib only; no slice, no package code).
include(joinpath(@__DIR__, "voracle.jl"))

S(f, t, c, r, re) = (f, t, Symbol[c...], Symbol[r...], re)
RE, SS = true, false

const CASES = [
    (name = "A re_turnover_beside_ss_route [S -> P]", subs = [:S], prods = [:P], ref = :S,
     faces = [], steps = [
        S("E", "E(S)", [:S], [], RE), S("E(S)", "E(P)", [], [], RE), S("E", "E(P)", [:P], [], RE),
        S("E(S)", "F(P)", [], [], SS), S("F(P)", "E", [], [:P], SS)]),
    (name = "B mut/expand_bibi/453/flip2 [A + B -> P + Q]", subs = [:A, :B], prods = [:P, :Q],
     ref = :A, faces = [], steps = [
        S("E", "E(A)", [:A], [], RE), S("E", "E(P)", [:P], [], RE), S("E", "E(Q)", [:Q], [], RE),
        S("E(A)", "E(P,Q)", [:B], [], RE), S("E(P)", "E(P,Q)", [:Q], [], RE),
        S("E(Q)", "E(P,Q)", [:P], [], SS)]),
    (name = "C hand LDH-like trap [A + B -> P + Q]", subs = [:A, :B], prods = [:P, :Q], ref = :A,
     faces = [[:P, :Q], [:Q], [:P]], steps = [
        S("E", "E(A)", [:A], [], RE), S("E(A)", "E(A,B)", [:B], [], SS),
        S("E(A,B)", "E(P,Q)", [], [], SS), S("E(P)", "E(P,Q)", [:Q], [], RE),
        S("E(Q)", "E(P,Q)", [:P], [], RE), S("E", "E(P)", [:P], [], SS), S("E", "E(Q)", [:Q], [], RE)]),
    (name = "D enumerator flip child 1 (split E/E(B)) [A + B -> P + Q]", subs = [:A, :B],
     prods = [:P, :Q], ref = :A, faces = [[:P, :Q], [:Q], [:P]], steps = [
        S("E", "E(A)", [:A], [], SS), S("E", "E(B)", [:B], [], SS), S("E", "E(Q)", [:Q], [], RE),
        S("E(A)", "E(P,Q)", [:B], [], SS), S("E(B)", "E(B,Q)", [:Q], [], RE),
        S("E(B)", "E(P,Q)", [:A], [], SS), S("E(Q)", "E(B,Q)", [:B], [], RE),
        S("E(Q)", "E(P,Q)", [:P], [], SS)]),
    (name = "E spec 40 LDH, active conformation as a plain network", subs = [:NADH, :Pyr],
     prods = [:Lac, :NAD], ref = :NADH, faces = [[:Lac, :NAD], [:NAD], [:Lac]], steps = [
        S("E", "E(Lac)", [:Lac], [], SS), S("E", "E(NAD)", [:NAD], [], RE),
        S("E", "E(NADH)", [:NADH], [], SS), S("E", "E(Pyr)", [:Pyr], [], RE),
        S("E(Lac)", "E(Lac,NAD)", [:NAD], [], RE), S("E(Lac)", "E(Lac,NADH)", [:NADH], [], SS),
        S("E(NAD)", "E(Lac,NAD)", [:Lac], [], RE), S("E(NAD)", "E(NAD,Pyr)", [:Pyr], [], RE),
        S("E(NADH)", "E(Lac,NADH)", [:Lac], [], SS), S("E(NADH)", "E(NADH,Pyr)", [:Pyr], [], RE),
        S("E(Pyr)", "E(NAD,Pyr)", [:NAD], [], RE), S("E(Pyr)", "E(NADH,Pyr)", [:NADH], [], SS),
        S("E(Lac,NAD)", "E(NADH,Pyr)", [], [], SS)]),
    (name = "F spec 14 [S1 + S2 -> P] (reverse face)", subs = [:S1, :S2], prods = [:P], ref = :S1,
     faces = [[:S1, :S2]], steps = [
        S("E", "E(P)", [:P], [], RE), S("E", "E(S1)", [:S1], [], SS), S("E", "E(S2)", [:S2], [], RE),
        S("E(P)", "E(S1,S2)", [], [], SS), S("E(S1)", "E(S1,S2)", [:S2], [], RE),
        S("E(S2)", "E(S1,S2)", [:S1], [], RE)]),
    (name = "G hand ss_dead_end_at_face [S -> P + Q]", subs = [:S], prods = [:P, :Q], ref = :S,
     faces = [[:P, :Q]], steps = [
        S("E", "E(S)", [:S], [], RE), S("E(S)", "E(P,Q)", [], [], SS),
        S("E(Q)", "E(P,Q)", [:P], [], RE), S("E", "E(Q)", [:Q], [], SS),
        S("E(P)", "E(P,Q)", [:Q], [], SS)]),
    (name = "H hand no_turnover_path [S -> P]", subs = [:S], prods = [:P], ref = :S,
     faces = [], steps = [
        S("E", "E(S)", [:S], [], RE), S("E(S)", "E(P)", [], [], SS), S("E(P)", "F", [], [:P], SS)]),
    (name = "I enumerator flip child 9 (trap E(Q)) [A + B -> P + Q]", subs = [:A, :B],
     prods = [:P, :Q], ref = :A, faces = [[:P, :Q]], steps = [
        S("E", "E(A)", [:A], [], SS), S("E", "E(P)", [:P], [], RE), S("E", "E(Q)", [:Q], [], SS),
        S("E(A)", "E(A,B)", [:B], [], RE), S("E(A)", "E(A,Q)", [:Q], [], SS),
        S("E(P)", "E(B,P)", [:B], [], RE), S("E(P)", "E(P,Q)", [:Q], [], RE),
        S("E(Q)", "E(A,Q)", [:A], [], SS), S("E(Q)", "E(P,Q)", [:P], [], RE),
        S("E(A,B)", "E(P,Q)", [], [], SS)]),
]

relspread(v) = (lo = minimum(v); hi = maximum(v); (hi - lo) / max(abs(hi), abs(lo), big(1e-300)))

function run_case(cs, rng; npar = 5, verbose_par = 1)
    println("\n=== ", cs.name)
    foreach(s -> println("    ", s[1], join([" + $x" for x in s[3]]), s[5] ? " ⇌ " : " <--> ",
                         s[2], join([" + $x" for x in s[4]])), cs.steps)
    mets = unique!(vcat(cs.subs, cs.prods, [x for s in cs.steps for x in vcat(s[3], s[4])]))
    for ip in 1:npar
        steps, mu = consistent(cs.steps, rng; mets)
        base = Dict(x => lu10(rng, -1, 1) for x in mets)
        sp1 = ones(BigFloat, length(steps)); sp2 = [lu10(rng, -3, 3) for _ in steps]
        vp = [rate_at(steps, base, big(10)^e, sp1, cs.ref)[1] for e in (10, 20, 40)]
        vp2 = rate_at(steps, base, big(10)^40, sp2, cs.ref)[1]
        ip <= verbose_par && @printf("  param %d, positive point: Λ=1e10,1e20,1e40: %s; Λ=1e40 other speeds: %s\n",
                                     ip, join(fmt.(vp), ", "), fmt(vp2))
        ip > verbose_par && @printf("  param %d, positive point Λ=1e40 ratio to Λ=1e20: %s\n", ip,
                                    fmt(vp[3] / vp[2]))
        for zs in cs.faces
            il = interior_limits(steps, base, zs, cs.ref, rng)
            fl = face_limits(steps, base, zs, cs.ref, rng)
            lims = [r[end] for r in il]
            law = lims[1]
            dir = relspread(lims)
            facev = fl.vals
            ok = all(!isnan, facev)
            spread = ok ? relspread(facev) : big(NaN)
            dev = ok ? maximum(abs.(facev .- law)) / max(abs(law), big(1e-300)) : big(NaN)
            if ip <= verbose_par
                @printf("    face %s = 0: interior limit (3 directions, ε=1e-20→1e-30): %s\n", join(zs, ","),
                        join(["$(fmt(r[1]))→$(fmt(r[2]))" for r in il], " | "))
                @printf("      exact face, 20 speed patterns in 10^[-3,3]: min %s max %s (terminal classes %d, Λ conv %s)\n",
                        fmt(minimum(facev)), fmt(maximum(facev)), fl.nterm, fmt(fl.conv))
            end
            @printf("    [param %d] face %-8s direction spread %-10s face-speed spread %-10s max |face-law|/|law| %s\n",
                    ip, join(zs, ","), fmt(dir), fmt(spread), fmt(dev))
        end
    end
end

rng = Xoshiro(20261011)
for cs in CASES
    run_case(cs, rng)
end
