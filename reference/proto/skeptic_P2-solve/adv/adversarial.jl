# ABOUTME: Skeptic of probe P2: runs Thermo.solve on adversarial mechanisms (package-built and
# ABOUTME: hand-built) and checks counts, semantics, thermodynamics, flux-free flags, rank, scrambles.
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/skeptic_P2-solve/adv/AdvLib.jl")
using Random, Printf, Serialization, LinearAlgebra
const K = joinpath(PROTO, "skeptic_P2-solve")
setprecision(BigFloat, 512)
rel(a, b) = abs(a - b) / max(abs(b), big(10)^-300)

# ─── hand-built mechanisms the package rejects ─────────────────────────────
M(n, r) = Pop.PMet(n, r)
const S_ = M(:S, :S); const P_ = M(:P, :P); const NONE = Pop.PMet[]
sp(conf, b...) = Pop.PSpecies(collect(b), conf, Symbol[], Symbol[],
                              Symbol(conf, join(string.(getfield.(collect(b), :name)))))
mk(id, species, steps, groups; inh = Symbol[]) = Pop.PMech(; id = id, source = id,
    substrates = [:S], products = [:P], inhibitors = inh, metabolites = [:S; :P; inh],
    species, steps, groups)
hand = [
    # E(S) <--> E(P) twice, in two groups (the package: "a reaction belongs to one group")
    mk("hand:parallel_identical_steps", [sp(:E), sp(:E, S_), sp(:E, P_)],
       [Pop.PStep(1, 2, [S_], NONE, true), Pop.PStep(2, 3, NONE, NONE, false),
        Pop.PStep(2, 3, NONE, NONE, false), Pop.PStep(1, 3, [P_], NONE, true)],
       [[1], [2], [3], [4]]),
    # (E <--> F, F <--> G) in one group (the package: different conformational changes)
    mk("hand:group_different_conformational_changes",
       [sp(:E), sp(:E, S_), sp(:E, P_), sp(:F), sp(:G)],
       [Pop.PStep(1, 2, [S_], NONE, false), Pop.PStep(2, 3, NONE, NONE, false),
        Pop.PStep(1, 3, [P_], NONE, false), Pop.PStep(1, 4, NONE, NONE, false),
        Pop.PStep(4, 5, NONE, NONE, false), Pop.PStep(5, 1, NONE, NONE, false)],
       [[1], [2], [3], [4, 5], [6]]),
]
adv = deserialize(joinpath(K, "data", "adv_pkg.jls"))
mechs = [adv; hand]

"""FD rank of the PACKAGE's parameterization (fitted params, its dependents, Pass-2 rename)."""
function fd_rank_pkg(m, rng; h = big(10)^-28)
    fitted = m.fitted_params; p = length(fitted)
    npts = 3p + 8
    vals0 = Dict{Symbol, BigFloat}(x => Oracle.logu(rng) for x in fitted)
    vals0[:Keq] = Oracle.logu(rng)
    cs = [Dict{Symbol, BigFloat}(x => Oracle.logu(rng) for x in m.metabolites) for _ in 1:npts]
    function consts(vals)
        v = copy(vals)
        for (d, mono) in m.dependent_monomials
            v[d] = prod((v[x]^BigFloat(e) for (x, e) in mono); init = big(1))
        end
        Oracle.constants_from_names(m, v; rename = Dict(m.pass2_rename))
    end
    J = zeros(BigFloat, npts, p)
    for (col, x) in enumerate(fitted), sg in (1, -1)
        v = copy(vals0); v[x] *= exp(sg * h)
        k = consts(v)
        for r in 1:npts; J[r, col] += sg * Oracle.steady_state_v(m, k, cs[r]) / (2h); end
    end
    k0 = consts(vals0)
    for r in 1:npts; J[r, :] ./= abs(Oracle.steady_state_v(m, k0, cs[r])); end
    d = sort!(abs.(diag(qr(J, ColumnNorm()).R)); rev = true); d ./= d[1]
    count(>(big(10)^-20), d)
end

rng = Xoshiro(2026)
for m in mechs
    println("\n== ", m.source, "  (", m.id, ")")
    for (j, s) in enumerate(m.steps)
        g = findfirst(gr -> j in gr, m.groups)
        println("   ", j, " [g", g, "] ", Pop.step_text(m, s))
    end
    sol = try
        Thermo.solve(m)
    catch e
        println("   SOLVE ERROR: ", first(sprint(showerror, e), 300)); continue
    end
    np = length(sol.indep)
    a = Adv.true_count(m; semantics = :package)
    b = Adv.true_count(m; semantics = :probe)
    allp = length(Thermo.solve(m; mode = :all).indep)
    pkg = isempty(m.samples) ? "-" : string(length(m.fitted_params))
    @printf("   counts: probe %d | package fitted %s | potentials, package groups %d%s | potentials, probe groups %d%s | AllParameters probe %d, classes %d\n",
            np, pkg, a.count, a.contradiction ? " CONTRADICTION" : "", b.count,
            b.contradiction ? " CONTRADICTION" : "", allp, a.allp)
    for d in sol.deps
        d.label === :kinetic_group && d.fold && continue
        println("   dep ", Thermo.colname(sol, d.col), " = Keq^", d.keq, " * ",
                join([string(Thermo.colname(sol, f), "^", e) for (f, e) in d.terms], " * "),
                "  [", d.label, d.fold ? ", fold" : "", "]")
    end
    es = Adv.exponent_stats(sol)
    println("   exponents: non-integer ", es.nonint, ", max |e| ", es.maxexp, ", max RREF entry ", es.maxint)
    # the probe's constants against the package's stored-orientation ties, and detailed balance
    k, keq, _ = Adv.random_constants(sol, rng)
    tie = Adv.ties_hold(m, k)
    v0 = Oracle.steady_state_v(m, k, Oracle.concs_at_ratio(m, rng, keq))
    v1 = Oracle.steady_state_v(m, k, Oracle.concs_at_ratio(m, rng, keq / 10))
    @printf("   probe constants: package ties hold %s; detailed balance |v(Keq)|/|v(Keq/10)| = %.1e\n",
            tie, Float64(abs(v0) / abs(v1)))
    # against the package's own constants
    if !isempty(m.samples)
        wb = big(0); wv = big(0)
        for s in m.samples
            vals = Dict{Symbol, BigFloat}(p => BigFloat(v) for (p, v) in s.params)
            for (d, mono) in m.dependent_monomials
                vals[d] = prod((vals[x]^BigFloat(e) for (x, e) in mono); init = big(1))
            end
            kp = Oracle.constants_from_names(m, vals; rename = Dict(m.pass2_rename))
            bare = Thermo.bare_values(sol, kp)
            x = Thermo.all_values(sol, Dict(i => bare[i] for i in sol.indep), vals[:Keq])
            wb = max(wb, maximum(rel(x[i], bare[i]) for i in keys(bare)))
            c = Dict{Symbol, BigFloat}(q => BigFloat(v) for (q, v) in s.concs)
            vp = Oracle.steady_state_v(m, kp, c)
            vm = Oracle.steady_state_v(m, Thermo.stored_constants(sol, x), c)
            wv = max(wv, rel(vm, BigFloat(s.v)), rel(vp, BigFloat(s.v)))
        end
        @printf("   vs package: worst bare-constant rebuild %.1e; worst oracle v (package constants or probe rebuild) vs package v %.1e\n",
                Float64(wb), Float64(wv))
    end
    # flux-free flags against numeric step fluxes at Γ = 10·Keq (probe constants, groups kept)
    z = Thermo.zero_flux_ss(m)
    br = Adv.plain_bridges(m)
    J, v = Adv.step_fluxes(m, k, Oracle.concs_at_ratio(m, rng, 10keq))
    sc = maximum(abs(x) for x in J if x !== nothing; init = abs(v))
    fl = [x === nothing ? "RE" : @sprintf("%.0e", Float64(abs(x) / sc)) for x in J]
    println("   flux-free flags ", findall(z), " | plain bridges ", findall(br .& .!getfield.(m.steps, :re)),
            " | |J|/max per step ", join(fl, " "))
    miss = [j for (j, x) in enumerate(J) if x !== nothing && !z[j] && abs(x) < big(10)^-30 * sc]
    falsepos = [j for (j, x) in enumerate(J) if x !== nothing && z[j] && abs(x) > big(10)^-30 * sc]
    isempty(miss) || println("   ZERO FLUX BUT NOT FLAGGED: steps ", miss)
    isempty(falsepos) || println("   FLAGGED BUT CARRIES FLUX: steps ", falsepos)
    # identifiable rank
    r = Adv.fd_rank(m, sol, rng)
    line = @sprintf("   fd rank (probe parameterization) %d (kept %.1e, dropped %.1e) vs probe count %d", r[1], r[2], r[3], np)
    if !isempty(m.samples)
        line *= @sprintf(" | fd rank (package parameterization) %d vs package count %d",
                         fd_rank_pkg(m, rng), length(m.fitted_params))
    end
    println(line)
    m2, nconv = Thermo.with_flux_free_as_re(m)
    if nconv > 0
        s2 = Thermo.solve(m2)
        r2 = Adv.fd_rank(m2, s2, rng)
        @printf("   flux-free groups made RE: %d -> count %d, fd rank %d\n", nconv, length(s2.indep), r2[1])
    end
    # scrambles
    sig = Adv.signature(sol); nf = 0
    for _ in 1:30
        s3 = try
            Thermo.solve(Adv.scramble(m, rng))
        catch e
            nf += 1; continue
        end
        Adv.signature(s3) == sig || (nf += 1)
    end
    println("   scrambles: ", nf, " of 30 differ")
end
