# ABOUTME: Skeptic of probe P2: follow-ups on the adversarial findings (Λ check of enumerated
# ABOUTME: phantoms, AllParameters vs the flux-free conversion, iso-group semantics, Keq^(1/2)).
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/skeptic_P2-solve/adv/AdvLib.jl")
using Random, Printf, Serialization
const K = joinpath(PROTO, "skeptic_P2-solve")
setprecision(BigFloat, 512)
rel(a, b) = abs(a - b) / max(abs(b), big(10)^-300)
rng = Xoshiro(777)
pops = Dict(p => Pop.load_population(p) for p in Pop.populations())
adv = Dict(m.source => m for m in deserialize(joinpath(K, "data", "adv_pkg.jls")))

println("— 1. enumerated count > rank: FD rank at Λ = 1e40 and 1e70, two random points each —")
for id in ("init_bibi_pp/61", "init_bibi_pp/62", "expand_bibi/303", "expand_bibi/304",
           "expand_bibi/305", "handwritten/23")
    p, i = split(id, "/"); m = pops[p][parse(Int, i)]
    sol = Thermo.solve(m)
    r1 = Adv.fd_rank(m, sol, rng); r2 = Adv.fd_rank(m, sol, rng; Λ = big(10)^70)
    @printf("%-16s count %d (package fitted %d)  rank Λ=1e40: %d (dropped %.0e)  Λ=1e70: %d (dropped %.0e)  %s\n",
            id, length(sol.indep), length(m.fitted_params), r1[1], r1[3], r2[1], r2[3], m.source)
end

println("\n— 2. flux-free conversion under AllParameters semantics (independent per-class constants) —")
for id in ("exchange_bridge", "dead_end_ss_binding", "exchange_in_wegscheider_cycle",
           "ss_in_re_segment", "inh_copy_iso_cycle", "dead_end_ss_triangle")
    m = adv[id]
    m2, nconv = Thermo.with_flux_free_as_re(m)
    worst = big(0); worst_c = big(0)
    for t in 1:3
        k, vals = Oracle.group_constants(m, rng)                 # AllParameters: no thermodynamics
        k2 = [m2.steps[j].re && !m.steps[j].re ?
              [Pop.is_binding(m.steps[j]) ? k[j][2] / k[j][1] : k[j][1] / k[j][2]] : k[j]
              for j in eachindex(k)]
        c = Dict{Symbol, BigFloat}(x => Oracle.logu(rng) for x in m.metabolites)
        va = Oracle.steady_state_v(m, k, c); vb = Oracle.steady_state_v(m2, k2, c)
        worst = max(worst, rel(vb, va))
        # the same with thermodynamically consistent constants (IndependentParameters)
        sol = Thermo.solve(m)
        kc, _, _ = Adv.random_constants(sol, rng)
        kc2 = [m2.steps[j].re && !m.steps[j].re ?
               [Pop.is_binding(m.steps[j]) ? kc[j][2] / kc[j][1] : kc[j][1] / kc[j][2]] : kc[j]
               for j in eachindex(kc)]
        worst_c = max(worst_c, rel(Oracle.steady_state_v(m2, kc2, c), Oracle.steady_state_v(m, kc, c)))
    end
    z = Thermo.zero_flux_ss(m); b = Adv.plain_bridges(m)
    @printf("%-32s flux-free steps %s (plain bridges %s): v(RE-converted)/v(SS) rel diff, consistent constants %.0e, independent per-class constants %.0e\n",
            id, findall(z), findall(z .& b), Float64(worst_c), Float64(worst))
end

println("\n— 3. iso groups: who reproduces the package's v? —")
for id in ("iso_group_opposite_chem", "iso_group_mixed_tier")
    m = adv[id]; sol = Thermo.solve(m)
    for s in m.samples[1:2]
        vals = Dict{Symbol, BigFloat}(p => BigFloat(v) for (p, v) in s.params)
        for (d, mono) in m.dependent_monomials
            vals[d] = prod((vals[x]^BigFloat(e) for (x, e) in mono); init = big(1))
        end
        kp = Oracle.constants_from_names(m, vals; rename = Dict(m.pass2_rename))
        c = Dict{Symbol, BigFloat}(q => BigFloat(v) for (q, v) in s.concs)
        bare = Thermo.bare_values(sol, kp)
        x = Thermo.all_values(sol, Dict(i => bare[i] for i in sol.indep), vals[:Keq])
        @printf("%-26s oracle(package constants) vs package v %.0e | oracle(probe rebuild from the same independents) vs package v %.2f\n",
                id, Float64(rel(Oracle.steady_state_v(m, kp, c), BigFloat(s.v))),
                Float64(rel(Oracle.steady_state_v(m, Thermo.stored_constants(sol, x), c), BigFloat(s.v))))
    end
    cs = Thermo.canonical_steps(m); al = Thermo.group_alignment(m, cs, m.groups)
    for g in m.groups
        length(g) < 2 && continue
        for j in g
            l, f = al[j]
            println("   member ", Pop.step_text(m, m.steps[j]), " (canonical ", cs[j].key, ") tied to lead ",
                    cs[l].key, f ? " REVERSED" : " as is", " | package: stored as written above, tied as stored")
        end
    end
end

println("\n— 4. hand-built: Keq^(1/2), an RE dead-end chain, and the fuzz dump —")
M(n, r) = Pop.PMet(n, r)
S_ = M(:S, :S); P_ = M(:P, :P); I_ = M(:I, :Inh); NONE = Pop.PMet[]
sp(conf, b...) = Pop.PSpecies(collect(b), conf, Symbol[], Symbol[],
                              Symbol(conf, join(string.(getfield.(collect(b), :name)))))
mk(id, species, steps, groups; inh = Symbol[]) = Pop.PMech(; id, source = id, substrates = [:S],
    products = [:P], inhibitors = inh, metabolites = [:S; :P; inh], species, steps, groups)
# E(S) -x-> F(S) -fused release-> E + P -x-> F -bind S, F->E-> E(S): x (E->F) twice per turnover
half_keq = mk("hand:half_keq", [sp(:E), sp(:F), sp(:E, S_), sp(:F, S_)],
    [Pop.PStep(3, 4, NONE, NONE, false),           # E(S) <--> F(S)   (group x)
     Pop.PStep(1, 2, NONE, NONE, false),           # E <--> F         (group x)
     Pop.PStep(4, 1, NONE, [P_], false),           # F(S) <--> E + P  (fused release)
     Pop.PStep(2, 3, [S_], NONE, false)],          # F + S <--> E(S)  (binding, F -> E)
    [[1, 2], [3], [4]])
# RE uni-uni plus an RE dead-end chain E + I ⇌ E(I) ⇌ F(I)
re_chain = mk("hand:re_dead_end_chain", [sp(:E), sp(:E, S_), sp(:E, P_), sp(:E, I_), sp(:F, I_)],
    [Pop.PStep(1, 2, [S_], NONE, true), Pop.PStep(2, 3, NONE, NONE, false),
     Pop.PStep(1, 3, [P_], NONE, true), Pop.PStep(1, 4, [I_], NONE, true),
     Pop.PStep(4, 5, NONE, NONE, true)], [[1], [2], [3], [4], [5]]; inh = [:I])
fdump = isfile(joinpath(K, "data", "fuzz_dump.jls")) ? deserialize(joinpath(K, "data", "fuzz_dump.jls")) : Pop.PMech[]
for m in [half_keq; re_chain; fdump]
    println("== ", m.id)
    for (j, s) in enumerate(m.steps)
        println("   ", j, " [g", findfirst(g -> j in g, m.groups), "] ", Pop.step_text(m, s))
    end
    sol = Thermo.solve(m)
    tc = Adv.true_count(m)
    for d in sol.deps
        println("   dep ", Thermo.colname(sol, d.col), " = Keq^", d.keq, " * ",
                join([string(Thermo.colname(sol, f), "^", e) for (f, e) in d.terms], " * "),
                "  [", d.label, d.fold ? ", fold" : "", "]")
    end
    k, keq, _ = Adv.random_constants(sol, rng)
    v0 = Oracle.steady_state_v(m, k, Oracle.concs_at_ratio(m, rng, keq))
    v1 = Oracle.steady_state_v(m, k, Oracle.concs_at_ratio(m, rng, keq / 10))
    z = Thermo.zero_flux_ss(m)
    J, v = Adv.step_fluxes(m, k, Oracle.concs_at_ratio(m, rng, 10keq))
    sc = maximum(abs(x) for x in J if x !== nothing; init = abs(v))
    m2, nconv = Thermo.with_flux_free_as_re(m)
    s2 = nconv > 0 ? Thermo.solve(m2) : sol
    r = Adv.fd_rank(m2, s2, rng)
    @printf("   count probe %d, potentials %d; detailed balance %.0e; flux-free %s; |J|/max %s; after conversion count %d, FD rank %d (dropped %.0e)\n",
            length(sol.indep), tc.count, Float64(abs(v0) / abs(v1)), findall(z),
            join([x === nothing ? "RE" : @sprintf("%.0e", Float64(abs(x) / sc)) for x in J], " "),
            length(s2.indep), r[1], r[3])
end
