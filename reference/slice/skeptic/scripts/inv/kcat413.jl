# ABOUTME: Investigates the kcat difference on P2's fuzz/413: structure, cores, the slice's kcat
# ABOUTME: and its own law at growing substrates, and the oracle on both cores.
include(joinpath(@__DIR__, "..", "spotlib.jl"))
r = only(x for x in skeptic_records("p2_fuzz_dump_flip") if x.p.id == "fuzz/413")
p, c, m = r.p, r.c, r.c.mech
fwd = S.orient(m)
for k in eachindex(m.edges)
    println("  edge $k g$(m.group[k])  ", S.step_text(m, S.oriented(m, k, fwd)))
end
core = my_core(p)
println("slice core ", Int.(S.catalytic_core(m)), "; my core ",
        [Int(core[findfirst(==(k), c.edge_of)]) for k in eachindex(m.edges)])
rng = Xoshiro(4)
d = S.derive(m)
sys = my_system(p, Symbol[]); mine = my_solution(sys, rng)
keq = big(3)
vals = my_sample(sys, mine, rng, log(keq))
kA = step_values(p, sys, vals, :A)
params = slice_params(c, d.solution, kA, kA, keq, rng)
println("slice kcat ", Float64(S.kcat(m, params)))
em = S.compile(m)
for X in (big(10)^5, big(10)^10, big(10)^20, big(10)^40)
    cc = Dict{Symbol, BigFloat}(x => x in m.substrates ? X : BigFloat(10)^-30 for x in S.metabolites(m))
    vl = S.rate_equation(em, (; (x => cc[x] for x in S.metabolites(m))...), params)
    vfull = Oracle.steady_state_v(p, kA, cc; Λ = big(10)^100)
    _, vcore, _ = my_network(p, kA, core, sort!(unique([f for (j, s) in enumerate(p.steps)
                               if core[j] for f in (s.from, s.to)])), cc; Λ = big(10)^100)
    @printf("X=1e%d: law %.10e  oracle full %.10e  oracle my core %.10e\n", round(Int, log10(X)),
            vl, vfull, vcore)
end
