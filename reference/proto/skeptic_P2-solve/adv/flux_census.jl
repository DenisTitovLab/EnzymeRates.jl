# ABOUTME: Skeptic of probe P2: numeric check of Thermo.zero_flux_ss on every plain mechanism and
# ABOUTME: on the second adversarial batch: the flagged steps' fluxes, unflagged zero fluxes.
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/skeptic_P2-solve/adv/AdvLib.jl")
using Random, Printf, Serialization
setprecision(BigFloat, 512)
const K = joinpath(PROTO, "skeptic_P2-solve")
rng = Xoshiro(8080)
nm = 0; nflag = 0; fp = String[]; miss = String[]; worst_flag = big(0)
for pop in Pop.populations(), m in Pop.load_population(pop)
    m.allosteric && continue
    global nm += 1
    sol = Thermo.solve(m); z = Thermo.zero_flux_ss(m)
    global nflag += count(z)
    for t in 1:2
        k, keq, _ = Adv.random_constants(sol, rng)
        J, v = Adv.step_fluxes(m, k, Oracle.concs_at_ratio(m, rng, 10keq))
        sc = maximum(abs(x) for x in J if x !== nothing; init = abs(v))
        for (j, x) in enumerate(J)
            x === nothing && continue
            r = abs(x) / sc
            z[j] && (global worst_flag = max(worst_flag, r))
            z[j] && r > big(10)^-30 && push!(fp, "$(m.id) step $j $(Float64(r))")
            !z[j] && r < big(10)^-30 && push!(miss, "$(m.id) step $j $(Pop.step_text(m, m.steps[j]))")
        end
    end
end
@printf("plain mechanisms %d, flagged flux-free SS steps %d; worst |J|/max over flagged steps %.1e\n",
        nm, nflag, Float64(worst_flag))
println("flagged but carrying flux: ", length(fp)); foreach(println, fp[1:min(end, 10)])
println("unflagged with zero flux: ", length(unique(miss))); foreach(println, unique(miss)[1:min(end, 10)])

println("\n— second adversarial batch (package-built) —")
for m in deserialize(joinpath(K, "data", "adv_pkg2.jls"))
    println("== ", m.source)
    for (j, s) in enumerate(m.steps)
        println("   ", j, " [g", findfirst(g -> j in g, m.groups), "] ", Pop.step_text(m, s), "   ", m.step_names[j])
    end
    sol = Thermo.solve(m)
    println("   package fitted ", m.fitted_params)
    for d in m.dependents; println("   package dep ", d); end
    for d in sol.deps
        println("   probe dep ", Thermo.colname(sol, d.col), " = Keq^", d.keq, " * ",
                join([string(Thermo.colname(sol, f), "^", e) for (f, e) in d.terms], " * "),
                "  [", d.label, d.fold ? ", fold" : "", "]")
    end
    k, keq, _ = Adv.random_constants(sol, rng)
    J, v = Adv.step_fluxes(m, k, Oracle.concs_at_ratio(m, rng, 10keq))
    sc = maximum(abs(x) for x in J if x !== nothing; init = abs(v))
    wb = big(0)
    for s in m.samples
        vals = Dict{Symbol, BigFloat}(p => BigFloat(x) for (p, x) in s.params)
        for (d, mono) in m.dependent_monomials
            vals[d] = prod((vals[x]^BigFloat(e) for (x, e) in mono); init = big(1))
        end
        kp = Oracle.constants_from_names(m, vals; rename = Dict(m.pass2_rename))
        bare = Thermo.bare_values(sol, kp)
        x = Thermo.all_values(sol, Dict(i => bare[i] for i in sol.indep), vals[:Keq])
        wb = max(wb, maximum(abs(x[i] - bare[i]) / abs(bare[i]) for i in keys(bare)))
    end
    r = Adv.fd_rank(m, sol, rng)
    @printf("   counts probe %d package %d; FD rank %d; rebuild of package constants %.0e; flux-free %s; |J|/max %s\n",
            length(sol.indep), length(m.fitted_params), r[1], Float64(wb), findall(Thermo.zero_flux_ss(m)),
            join([x === nothing ? "RE" : @sprintf("%.0e", Float64(abs(x) / sc)) for x in J], " "))
end
