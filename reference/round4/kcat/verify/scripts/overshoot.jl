# ABOUTME: Rechecks the "rate exceeds kcat" claims with an independent search: the probe's flagged
# ABOUTME: laws, the classes it says never overshoot, and hand-built branches with no abortive form.
include(joinpath(@__DIR__, "vlib.jl"))

include(joinpath(@__DIR__, "overshoot_lib.jl"))

t0 = time()
rng = Xoshiro(99)

# 1. hand-built branches without abortive forms
hand = [
    ("mixed ping-pong and sequential", [:A, :B], [:P, :Q],
     ["E + A <=> E(A)", "E(A) + B <--> E(A, B)", "E(A, B) <--> E(P, Q)", "E(Q) + P <=> E(P, Q)",
      "E + Q <=> E(Q)", "E(A) <--> E(P; +A -P)", "E(; +A -P) + P <=> E(P; +A -P)",
      "E(; +A -P) + B <--> E(B; +A -P)", "E(B; +A -P) <--> E(Q)"]),
    ("isomer branch at E(A)", [:A, :B], [:P, :Q],
     ["E + A <=> E(A)", "E(A) + B <--> E(A, B)", "E(A, B) <--> E(P, Q)", "E(Q) + P <=> E(P, Q)",
      "E + Q <=> E(Q)", "E(A) <--> F(A)", "F(A) + B <--> F(A, B)", "F(A, B) <--> E(P, Q)"]),
    ("textbook random-order SS bi-bi", [:A, :B], [:P, :Q],
     ["E + A <--> E(A)", "E + B <--> E(B)", "E + Q <=> E(Q)", "E(A) + B <--> E(A, B)",
      "E(B) + A <--> E(A, B)", "E(Q) + P <=> E(P, Q)", "E(A, B) <--> E(P, Q)"]),
]
println("== hand-built mechanisms (300 parameter points, constants 1e-3..1e3)")
for (name, subs, prods, lines) in hand
    m = mech(subs, prods, lines)
    cl = overshoot_class(m)
    r = test_law(m, 300, rng; plo = -3.0, phi = 3.0)
    println(name, ": class (SS on non-turnover cycle, abortive on it, SS at abortive) = ", cl)
    foreach(l -> println("    ", l), steps_text(m))
    if r.ok
        @printf("    max sup/kcat %.4g (kcat %.4g, sup %.4g at log10 conc %s); oracle vs law at argmax %.1e\n",
                r.ratio, r.kc, r.sup, string(round.(r.at; digits = 2)), r.orc)
    else
        println("    not tested: ", r.why)
    end
    flush(stdout)
end

# 2. the probe's flagged laws (sup > 1.01 kcat) and the plain classes it says never overshoot
rowsP = [x for l in ("core", "expand_terter", "terter_sample")
         for x in deserialize(joinpath(KDATA, "rows_$l.jls")) if x.ok && x.cls === :finite]
flag = Set((x.set, x.index) for x in rowsP if x.ratio > 1.01)
println("\nprobe rows: finite ", length(rowsP), ", flagged ", length(flag))
recs = Dict{Tuple{String, Int}, Any}()
for set in CORE
    foreach(r -> recs[(r.set, r.index)] = r.m, core_records(set))
end
foreach(r -> recs[(r.set, r.index)] = r.m, pop_records("expand_terter"))
picks = deserialize(joinpath(KDATA, "terter_sample.jls"))
for set in unique(first.(picks))
    foreach(r -> recs[(r.set, r.index)] = r.m, pop_records(set, [i for (s, i, _) in picks if s == set]))
end
fl = collect(flag)
sample_flag = fl[randperm(rng, length(fl))[1:60]]
again = 0
worst_orc = 0.0
for k in sample_flag
    r = test_law(recs[k], 5, rng)
    r.ok || continue
    global again += r.ratio > 1.01
    isnan(r.orc) || (global worst_orc = max(worst_orc, r.orc))
end
@printf("flagged laws re-searched (60 random, 5 points each): overshoot again in %d; oracle vs law at argmax worst %.1e\n",
        again, worst_orc)
flush(stdout)

# plain finite laws by class, then a harder search on the classes the probe says never overshoot
plain = [(k, recs[k]) for k in ((x.set, x.index) for x in rowsP if !x.allo)]
cls = Dict(k => overshoot_class(m) for (k, m) in plain)
tab = Dict{Any, Vector{Int}}()
for (k, m) in plain
    v = get!(tab, cls[k], [0, 0]); v[1] += 1; v[2] += k in flag
end
println("plain finite laws by my class (SS on nt cycle, abortive on it, SS at abortive): [laws, probe flagged]")
foreach(p -> println("   ", p), sort!(collect(tab); by = string))
quiet = [(k, m) for (k, m) in plain if !cls[k][3] && !(k in flag)]
pick = quiet[randperm(rng, length(quiet))[1:min(500, end)]]
hits = Any[]
for (k, m) in pick
    r = test_law(m, 6, rng; plo = -3.0, phi = 3.0)
    r.ok && r.ratio > 1.01 && push!(hits, (k, cls[k], r))
end
@printf("classes without SS step at an abortive form: %d laws, %d searched (6 points, constants 1e-3..1e3, box 1e-8..1e8): sup > 1.01 kcat in %d\n",
        length(quiet), length(pick), length(hits))
for (k, c, r) in hits[1:min(4, end)]
    @printf("   %s/%d class %s ratio %.4g (kcat %.4g, sup %.4g), oracle %.1e\n", k[1], k[2], string(c),
            r.ratio, r.kc, r.sup, r.orc)
    foreach(l -> println("      ", l), steps_text(recs[k]))
end
# the flagged ones' classes under my definition
println("probe-flagged plain laws by my class: ",
        (d = Dict{Any, Int}(); foreach(k -> haskey(cls, k) && (d[cls[k]] = get(d, cls[k], 0) + 1), flag); d))
@printf("done in %.0f s, maxrss %.0f MB\n", time() - t0, Sys.maxrss() / 2^20)
