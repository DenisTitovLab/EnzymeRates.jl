# ABOUTME: specs/40 with the skeptic's consistent sample: law against the oracle from moderate to
# ABOUTME: saturating substrates, products at 0, and the slice's kcat at the same values.
include(joinpath(@__DIR__, "..", "spotlib.jl"))
r = only(x for x in skeptic_records("specs") if x.p.id == "specs/40")
p, c, m = r.p, r.c, r.c.mech
rng = Xoshiro(7)
d = S.derive(m)
tags = p.group_tags
sys = my_system(p, tags); mine = my_solution(sys, rng)
keq = big(3)
vals = my_sample(sys, mine, rng, log(keq))
kA, kI = step_values(p, sys, vals, :A), step_values(p, sys, vals, :I)
params = slice_params(c, d.solution, kA, kI, keq, rng)
flip, keepI = my_inactive(p, tags)
em = S.compile(m)
println("flip species ", p.species[flip].name, "; inactive steps ", findall(keepI))
for X in (big(10)^0, big(10)^5, big(10)^10, big(10)^20, big(10)^40), prods in (big(0), big(1))
    cc = Dict{Symbol, BigFloat}(x => x in m.substrates ? X : prods for x in S.metabolites(m))
    vo = my_mwc(p, kA, kI, keepI, flip, params.L, big(1), big(1), m.multiplicity, cc;
                Λ = big(10)^100)
    vl = S.rate_equation(em, (; (x => cc[x] for x in S.metabolites(m))...), params)
    @printf("X=1e%-3d products %d: oracle %.12e  law %.12e  rel %.1e\n",
            round(Int, log10(X)), Int(prods), vo, vl, relerr(vl, vo))
end
println("slice kcat ", Float64(S.kcat(m, params)))
println("-- products tending to 0 at X = 1e10")
for e in (0, 5, 10, 20, 40, -1)
    pr = e < 0 ? big(0) : BigFloat(10)^(-e)
    cc = Dict{Symbol, BigFloat}(x => x in m.substrates ? big(10)^10 : pr for x in S.metabolites(m))
    xA, JA, locA = my_network(p, kA, trues(length(p.steps)), collect(eachindex(p.species)), cc;
                              Λ = big(10)^100)
    formsI = sort!(unique([flip; [f for (j, s) in enumerate(p.steps) if keepI[j] for f in (s.from, s.to)]]))
    xI, JI, locI = my_network(p, kI, keepI, formsI, cc; Λ = big(10)^100)
    vl = S.rate_equation(em, (; (x => cc[x] for x in S.metabolites(m))...), params)
    @printf("products %s: πA %.6e JA %.6e πI %.6e JI %.6e | law %.10e\n",
            e < 0 ? "0" : "1e-$e", xA[locA[flip]], JA, xI[locI[flip]], JI, vl)
    if e < 0
        for (i, f) in enumerate(eachindex(p.species))
            @printf("   A %-22s %.3e\n", p.species[f].name, xA[i])
        end
    end
end
