# ABOUTME: R5.1 agreement and cost: the slice (U) and law-type variants A, B, C on the census sample:
# ABOUTME: equal laws and outputs, solve+law times interleaved per mechanism, live bytes per law.
include(joinpath(@__DIR__, "lib.jl"))
const VS = [(v, load_variant(v)) for v in ("U", "A", "B", "C")]
const SETS = split(get(ENV, "SETS", "spec,bibi,pingpong,twosites,uniuni"), ',')
const DOWN = get(ENV, "DOWN", "1") == "1"
const TAG = get(ENV, "TAG", join(SETS, "_"))
t_start = time()

rows = NamedTuple[]
nbuild = Dict{String, Int}()
for r in census_rows()
    group(r.set) in SETS || continue
    try
        push!(rows, (set = r.set, id = r.id, ms = [slice_mech(S, r.net) for (_, S) in VS]))
    catch e
        k = group(r.set) * ": " * first(sprint(showerror, e), 70)
        nbuild[k] = get(nbuild, k, 0) + 1
    end
end
println("julia ", VERSION, "; mechanisms: ", length(rows), " (", join(("$g $(count(r -> group(r.set) == g, rows))" for g in SETS), ", "),
        "); census rows the slice cannot build: ", isempty(nbuild) ? "none" : nbuild)
flush(stdout)

"""Everything a variant derives for `m`: the law, its printed string, key and code."""
function outputs(S, m)
    sol = S.solve(m)
    L = S.rate_law(m, sol)
    sol, L, (lawdata(L), S.rate_equation_string(m, L, sol), S.canonical(L).key,
             S.body(S.program(m, sol, L)))
end

# ── agreement and timing, variants interleaved per mechanism ───────────────────────────
res = Dict(v => Dict(k => Dict{String, Vector{Float64}}() for k in (:solve, :law, :sl, :down))
           for (v, _) in VS)
diffs, errs, nbig, nallo = Dict{String, Int}(), Dict{String, Int}(), 0, 0
push_t!(v, kind, g, x) = push!(get!(res[v][kind], g, Float64[]), x)
for (k, r) in enumerate(rows)
    g = group(r.set)
    out = Any[]
    for (i, (v, S)) in enumerate(VS)
        o = try
            outputs(S, r.ms[i])
        catch e
            sprint(showerror, e)
        end
        push!(out, o)
    end
    if out[1] isa String
        all(o -> o isa String, out) || (diffs[g * ": error in U only"] =
                                        get(diffs, g * ": error in U only", 0) + 1)
        key = g * ": " * first(out[1], 60)
        errs[key] = get(errs, key, 0) + 1
        continue
    end
    for i in 2:length(VS)
        out[i] isa String && (diffs["$g: $(VS[i][1]) errors"] = get(diffs, "$g: $(VS[i][1]) errors", 0) + 1; continue)
        for (j, part) in enumerate(("law", "string", "key", "code"))
            out[i][3][j] == out[1][3][j] ||
                (diffs["$g: $(VS[i][1]) $part"] = get(diffs, "$g: $(VS[i][1]) $part", 0) + 1)
        end
    end
    L1 = out[1][2]
    global nbig += (L1 isa VS[1][2].Law{BigInt} || L1 isa VS[1][2].AlloLaw{BigInt})
    global nallo += L1 isa VS[1][2].AlloLaw
    gg = L1 isa VS[1][2].AlloLaw ? g * "/allo" : g * "/plain"
    order = circshift(1:length(VS), k)                  # rotate which variant runs first
    for i in order
        v, S = VS[i]
        m, sol, L = r.ms[i], out[i][1], out[i][2]
        ts = best(() -> S.solve(m))
        tl = best(() -> S.rate_law(m, sol))
        for gx in (g, gg)
            push_t!(v, :solve, gx, ts)
            push_t!(v, :law, gx, tl)
            push_t!(v, :sl, gx, ts + tl)
        end
        DOWN && push_t!(v, :down, g,
                        best(() -> (S.rate_equation_string(m, L, sol); S.canonical(L);
                                    S.body(S.program(m, sol, L))); reps = 2))
    end
    k % 250 == 0 && (println("  ", k, " done, ", round(Int, time() - t_start), " s"); flush(stdout))
end
println("derivation errors (same in every variant unless listed under differences): ",
        isempty(errs) ? "none" : errs)
println("differences against U: ", isempty(diffs) ? "none" : diffs)
println("laws that needed BigInt in U: ", nbig, "; allosteric laws: ", nallo)
gs = sort!(unique(reduce(vcat, [collect(keys(res["U"][:sl]))])))
for kind in (:sl, :solve, :law, :down), gx in [gs; "all"]
    for (v, _) in VS
        haskey(res["U"][kind], gx) || gx == "all" || continue
        t = gx == "all" ? reduce(vcat, [res[v][kind][x] for x in SETS if haskey(res[v][kind], x)]; init = Float64[]) :
                          res[v][kind][gx]
        isempty(t) && continue
        tu = gx == "all" ? reduce(vcat, [res["U"][kind][x] for x in SETS if haskey(res["U"][kind], x)]; init = Float64[]) :
                           res["U"][kind][gx]
        rr = t ./ tu
        @printf("%-5s %-15s %s  %s µs | ratio to U: median %.3f p90 %.3f max %.3f\n",
                kind, gx, v, stats(t), q(rr, 0.5), q(rr, 0.9), maximum(rr))
    end
end
serialize(joinpath(R5, "data", "timing_$(TAG).jls"), res)
flush(stdout)

# ── live bytes: one pass per variant over the same mechanisms ──────────────────────────
if get(ENV, "LIVE", "1") == "1"
    for (label, f) in (("law", (S, m) -> S.rate_law(m, S.solve(m))),
                       ("law+print+key", (S, m) -> (sol = S.solve(m); L = S.rate_law(m, sol);
                                                    S.rate_equation_string(m, L, sol); S.canonical(L))))
        for i in (1, 4, 3, 2)                                  # U, C, B, then A's BigInts
            v, S = VS[i]
            foreach(r -> try f(S, r.ms[i]) catch end, rows[1:min(end, 50)])
            GC.gc(true); GC.gc(true)
            lb0, g0 = Base.gc_live_bytes(), Base.gc_num()
            t = @elapsed for r in rows
                try
                    f(S, r.ms[i])
                catch
                end
            end
            g1 = Base.gc_num()
            GC.gc(true); GC.gc(true)
            lb1 = Base.gc_live_bytes()
            @printf("live bytes %-14s %s: %d mechanisms in %.1f s; after a full collection %+.2f KiB per mechanism (%+.1f MiB in all); %d collections, %d full; allocated %.1f MiB\n",
                    label, v, length(rows), t, (lb1 - lb0) / 1024 / length(rows), (lb1 - lb0) / 2^20,
                    g1.pause - g0.pause, g1.full_sweep - g0.full_sweep, (Base.gc_total_bytes(g1) - Base.gc_total_bytes(g0)) / 2^20)
            flush(stdout)
        end
    end
end
println("finished in ", round(Int, time() - t_start), " s")
