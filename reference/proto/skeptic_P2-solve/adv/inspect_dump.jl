# ABOUTME: Skeptic of probe P2: prints the fuzz mechanisms that failed a check (group alignment,
# ABOUTME: per-step fluxes, flux-free flags, rank). Usage: julia inspect_dump.jl <dump name>.
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/skeptic_P2-solve/adv/AdvLib.jl")
using Random, Printf, Serialization
setprecision(BigFloat, 512)
rng = Xoshiro(5)
for m in deserialize(joinpath(PROTO, "skeptic_P2-solve", "data", ARGS[1] * ".jls"))
    println("== ", m.id)
    cs = Thermo.canonical_steps(m)
    al = Thermo.group_alignment(m, cs, m.groups)
    for (j, s) in enumerate(m.steps)
        g = findfirst(gr -> j in gr, m.groups)
        l, f = al[j]
        println("   ", j, " [g", g, "] ", rpad(Pop.step_text(m, s), 30), " canonical ", cs[j].key,
                length(m.groups[g]) > 1 ? string("  → tied to step ", l, f ? " REVERSED" : " as canonical") : "")
    end
    sol = Thermo.solve(m)
    k, keq, _ = Adv.random_constants(sol, rng)
    println("   package-semantics ties hold: ", Adv.ties_hold(m, k), "; counts probe ", length(sol.indep),
            ", potentials (package groups) ", Adv.true_count(m).count,
            ", potentials (probe groups) ", Adv.true_count(m; semantics = :probe).count)
    z = Thermo.zero_flux_ss(m)
    for t in 1:2
        J, v = Adv.step_fluxes(m, k, Oracle.concs_at_ratio(m, rng, 10keq))
        sc = maximum(abs(x) for x in J if x !== nothing; init = abs(v))
        println("   |J|/max: ", join([x === nothing ? "RE" : @sprintf("%.0e", Float64(abs(x) / sc)) for x in J], " "),
                "   flux-free flags ", findall(z))
        k, keq, _ = Adv.random_constants(sol, rng)
    end
    m2, nconv = Thermo.with_flux_free_as_re(m)
    s2 = nconv > 0 ? Thermo.solve(m2) : sol
    println("   count after conversion ", length(s2.indep), ", FD rank ", Adv.fd_rank(m2, s2, rng)[1])
end
