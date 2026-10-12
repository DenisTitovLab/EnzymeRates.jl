# ABOUTME: Investigates the kcat difference on specs/40: structure, the slice's candidates per
# ABOUTME: corner, and the oracle per corner on the slice's and the skeptic's kcat steps.
include(joinpath(@__DIR__, "..", "spotlib.jl"))
r = only(x for x in skeptic_records(ARGS[2]) if x.p.id == ARGS[1])
p, c, m = r.p, r.c, r.c.mech
println(p.source, " tags ", m.tags, " n ", m.multiplicity, " sites ", m.sites)
fwd = S.orient(m)
for k in eachindex(m.edges)
    println("  edge $k g$(m.group[k]) $(m.tags[m.group[k]])  ", S.step_text(m, S.oriented(m, k, fwd)))
end
println("catalytic_core ", Int.(S.catalytic_core(m)), " kcat_mask ", Int.(S.kcat_mask(m)))
core = my_core(p)
core .|= [s.re && Pop.is_free(p.species[s.from]) && Pop.is_free(p.species[s.to]) for s in p.steps]
println("my core per edge ", [Int(core[findfirst(==(k), c.edge_of)]) for k in eachindex(m.edges)])
rng = Xoshiro(7)
d = S.derive(m)
tags = p.group_tags
sys = my_system(p, tags); mine = my_solution(sys, rng)
lkeq = log(big(3)); keq = big(3)
vals = my_sample(sys, mine, rng, lkeq)
kA, kI = step_values(p, sys, vals, :A), step_values(p, sys, vals, :I)
params = slice_params(c, d.solution, kA, kI, keq, rng)
println("params ", params)
vars, terms = S.kcat_candidates(m)
x = BigFloat[vars.conc[i] ? big(1) : getfield(params, s) for (i, s) in enumerate(vars.names)]
println("slice candidates: ", [Float64(S.evaluate(n, x) / S.evaluate(dd, x)) for (n, dd) in terms])
flip, keepI = my_inactive(p, tags)
regvals = Dict(k => getfield(params, k) for k in keys(params))
function at(cc, use)
    Λ = big(10)^100
    formsA = sort!(unique([f for (j, s) in enumerate(p.steps) if use[j] for f in (s.from, s.to)]))
    xA, JA, locA = my_network(p, kA, use, formsA, cc; Λ)
    keepIc = keepI .& use
    formsI = sort!(unique([flip; [f for (j, s) in enumerate(p.steps) if keepIc[j] for f in (s.from, s.to)]]))
    xI, JI, locI = my_network(p, kI, keepIc, formsI, cc; Λ)
    n = m.multiplicity
    wA = my_W(m, regvals, cc, :A) / xA[locA[flip]]^n
    wI = params.L * my_W(m, regvals, cc, :I) / xI[locI[flip]]^n
    (wA * JA + wI * JI) / (wA + wI)
end
mods = [x for x in S.metabolites(m) if !S.is_reactant(m, x)]
println("modifiers ", mods)
for X in (big(10)^20, big(10)^30, big(10)^40)
    for sat in S.all_subsets(mods)
        cc = Dict{Symbol, BigFloat}(x => x in m.substrates ? X : big(0) for x in S.metabolites(m))
        foreach(y -> cc[y] = sqrt(X), sat)
        sl = at(cc, S.kcat_mask(m)[[c.edge_of[j] for j in eachindex(p.steps)]])
        my = at(cc, core)
        full = at(cc, trues(length(p.steps)))
        @printf("X=1e%d sat=%-12s slice-mask %.10e  my-core %.10e  full %.10e\n",
                round(Int, log10(X)), join(sat, ","), sl, my, full)
    end
end
println("slice kcat ", S.kcat(m, params))
