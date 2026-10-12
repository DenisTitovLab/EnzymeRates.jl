# ABOUTME: Tests the probe's speed-free (dead-end) predicate against the numeric test on enumerated mechanisms of
# ABOUTME: reactions with a dead-end inhibitor; reports mismatches and compares the exact block predicate.
include(joinpath(@__DIR__, "SkLib.jl"))
using Random, Printf, Serialization
setprecision(BigFloat, Oracle.PREC)
pop = deserialize(joinpath(@__DIR__, "data", length(ARGS) >= 1 ? ARGS[1] : "enum_inh.jls"))
stats = Dict{String, Vector{Int}}()     # label => [n, probe≠thermo, exact≠thermo, exact_bare≠bare, errors]
examples = String[]
oracle_err = Dict{String, Vector{Float64}}()
for (n, m) in enumerate(pop)
    label = split(m.id, "/")[1]
    s = get!(stats, label, zeros(Int, 6))
    s[1] += 1
    L = try
        Engine.bare_rate_law(m; check = true)
    catch e
        s[5] += 1; continue
    end
    rng = Random.Xoshiro(17 + n)
    kc, _, _ = Oracle.consistent_constants(m, rng; Keq = Oracle.logu(rng))
    kb = Sk.bare_constants(m, rng)
    c = Sk.random_concs(m, rng)
    th = Sk.speed_free_numeric(L, m, kc, c)
    push!(get!(oracle_err, label, Float64[]), Sk.relerr(Sk.lawv(L, Sk.vals_of(L, m, kb, c)), Oracle.steady_state_v(m, kb, c)))
    ba = Sk.speed_free_numeric(L, m, kb, c)
    pt, pb = Sk.speed_free_exact(m)
    pr = L.stats[:dead_end]
    s[2] += pr != th; s[3] += pt != th; s[4] += pb != ba; s[6] += !isempty(th)
    if pr != th && length(examples) < 12
        push!(examples, "$(m.id): probe $(pr) thermo $(th) exact $(pt) | " *
              join([Pop.step_text(m, x) * (j in th ? " [speed-free]" : "") for (j, x) in enumerate(m.steps)], "; "))
    end
end
for (label, s) in sort!(collect(stats))
    @printf("%-18s n=%4d  with speed-free steps %4d | probe predicate ≠ numeric(thermo) %4d | exact block predicate ≠ numeric(thermo) %d | exact bridge predicate ≠ numeric(bare) %d | engine errors %d\n",
            label, s[1], s[6], s[2], s[3], s[4], s[5])
end
for (label, e) in sort!(collect(oracle_err))
    @printf("%-18s Engine (Hill == edge sum checked) vs oracle at bare constants: max rel %.1e over %d mechanisms\n", label, maximum(e), length(e))
end
println("\nexamples where the probe's predicate misses or overcounts:")
foreach(e -> println("  ", e), examples)
