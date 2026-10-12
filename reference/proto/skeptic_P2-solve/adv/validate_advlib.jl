# ABOUTME: Skeptic of probe P2: checks the probe's counts with an independent method (form and
# ABOUTME: metabolite potentials, exact) on every plain mechanism, and sanity-checks fd_rank.
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/skeptic_P2-solve/adv/AdvLib.jl")
using Random, Printf

n = 0; bad_pkg = String[]; bad_probe = String[]; bad_allp = String[]; contra = 0
semdiff = String[]; tiefail = String[]
t0 = time()
for pop in Pop.populations(), m in Pop.load_population(pop)
    m.allosteric && continue
    global n += 1
    sol = Thermo.solve(m)
    a = Adv.true_count(m; semantics = :package)
    b = Adv.true_count(m; semantics = :probe)
    np = length(sol.indep)
    a.count == length(m.fitted_params) || push!(bad_pkg, "$(m.id) potentials $(a.count) pkg $(length(m.fitted_params))")
    b.count == np || push!(bad_probe, "$(m.id) potentials $(b.count) probe $np")
    a.count == b.count || push!(semdiff, "$(m.id) package-semantics $(a.count) probe-semantics $(b.count)")
    a.allp == length(Thermo.solve(m; mode = :all).indep) || push!(bad_allp, m.id)
    (a.contradiction || b.contradiction) && (global contra += 1)
    k, _, _ = Adv.random_constants(sol, Xoshiro(n))
    Adv.ties_hold(m, k) || push!(tiefail, m.id)
end
@printf("plain mechanisms %d (%.1f s)\n", n, time() - t0)
println("potential-method count != package fitted count: ", length(bad_pkg)); foreach(println, bad_pkg)
println("potential-method count (probe's group semantics) != probe count: ", length(bad_probe)); foreach(println, bad_probe)
println("package vs probe group semantics give different counts: ", length(semdiff)); foreach(println, semdiff)
println("AllParameters count != class count: ", length(bad_allp))
println("contradictions: ", contra)
println("probe constants break the package's stored-orientation ties: ", length(tiefail)); foreach(println, tiefail)

# fd_rank sanity on mechanisms whose rank the probe reported
rng = Xoshiro(99)
pops = Dict(p => Pop.load_population(p) for p in ("specs", "handwritten", "init_bibi"))
for (p, i) in (("specs", 2), ("specs", 3), ("specs", 14), ("specs", 1), ("handwritten", 17),
               ("handwritten", 24), ("init_bibi", 56))
    m = pops[p][i]
    sol = Thermo.solve(m)
    t = @elapsed r = Adv.fd_rank(m, sol, rng)
    @printf("%-16s count %2d  fd_rank %2d  kept %.1e dropped %.1e  (%.1f s)\n", m.id,
            length(sol.indep), r[1], r[2], r[3], t)
end
